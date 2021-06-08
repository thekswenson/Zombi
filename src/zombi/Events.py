from .Interval import Interval

# Types:
T_EVENT = str
T_DIRECTION = bool
RIGHT = 0
LEFT = 1

# Event types:
DUP = 'D'
FER = 'T'
LOSS = 'L'
INV = 'I'
POS = 'P'
ORIG = 'O'

class GenomeEvent:
    """
    An rearrangement event. Meant to be used as a base class.

    Attributes
    ----------
    etype: str
        the type of event from {DUP, FER, LOSS, INV, POS, ORIG}
    lineage: str
        the lineage on which the event happened (pendant node name)
    time: float
        the time at which it happened
    """
    def __init__(self, etype: T_EVENT, lineage: str, time: float):
        self.etype: T_EVENT = etype
        self.lineage: str = lineage
        self.time: float = time

class EventTwoCuts(GenomeEvent):
    """
    An event with two cuts. Meant to be used as a base class.
    """
    def __init__(self, int1: Interval, int2: Interval, sc1: int, sc2: int,
                 direction: T_DIRECTION, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.before1: Interval = int1
        self.before2: Interval = int2
        self.sc1: int = sc1
        self.sc2: int = sc2
        self.direction: T_DIRECTION = direction

class Origination(GenomeEvent):
    pass

class Inversion(EventTwoCuts):
    """
    An Inversion event.
    """
    def __init__(self, *args, **kwargs):
        """
        Create an Inversion event.  When instantiating this you must also
        provide the arguments for EventTwoCuts and GenomeEvent.

        Parameters
        ----------
        int1 : Interval
            the old intergene at the left
        int2 : Interval
            the old intergene at the right
        sc1 : int
            the specific coordinate where the first interval is to be cut
        sc2 : int
            the specific coordinate where the second interval is to be cut
        """
        super().__init__(*args, **kwargs)
        self.after1: Interval = None
        self.after2: Interval = None

        self.setAfter()

    def setAfter(self):
        """
        Set the two intergenic regions that exist after the inversion.
        Consider intergenic regions I = `before1` and J = `before2` on either
        side of segement S:

            I S J

        I is split at `c1` into I1 I2 and J is split at `c2` into J1 J2. Then
        we get

            I1 I2 S J1 J2

        and the inversion produces

            I1 -J1 -S -I2 J2
        """
        lenI1 = self.sc1 - self.before1.sc1
        lenI2 = self.before1.sc2 - self.sc1
        lenJ1 = self.sc2 - self.before2.sc1
        lenJ2 = self.before2.sc2 - self.sc2

        sleftstart = self.before1.sc1
        sleftend = sleftstart + lenI1 + lenJ1
        tleftstart = self.before1.tc1
        tleftend = tleftstart + lenI1 + lenJ1

        srightstart = self.before2.sc1
        srightend = srightstart + lenI2 + lenJ2
        trightstart = self.before2.tc1
        trightend = trightstart + lenI2 + lenJ2

        position = self.before1.position
        self.after1 = Interval(tleftstart, tleftend,
                               sleftstart, sleftend, position, 'I')
        self.after2 = Interval(trightstart, trightend,
                               srightstart, srightend, position+1, 'I')

class TandemDup(EventTwoCuts):
    """
    A tandem duplication event.
    """
    def __init__(self, *args, **kwargs):
        """
        Create a TandemDuplication event.  When instantiating this you must also
        provide the arguments for EventTwoCuts and GenomeEvent.

        Parameters
        ----------
        int1 : Interval
            the old intergene at the left
        int2 : Interval
            the old intergene at the right
        """
        super().__init__(*args, **kwargs)
        self.after1: Interval = None
        self.after2: Interval = None
        self.after3: Interval = None

        self.setAfter()


    def setAfter(self):
        """
        Set the three intergenic regions that exist after the tandem
        duplication.  Consider intergenic regions I = `before1` and
        J = `before2` on either side of segement S:

            I S J

        I is split at `c1` into I0 I1 and J is split at `c2` into J0 J1. Then
        we get

            I0 I1 S J0 J1

        and the tandem duplication produces

            I0 I1 S J0 I1 S J0 J1
        """
        self.after1 = self.before1

        lenI2 = self.before1.sc2 - self.sc1
        lenS = self.before2.sc1 - self.before1.sc2
        lenJ1 = self.sc2 - self.before2.sc1
        lenJ2 = self.before2.sc2 - self.sc2
        if self.direction == RIGHT:
            scenterstart = self.before2.sc1
            scenterend = self.before2.sc1 + lenJ1 + lenI2
            tcenterstart = self.before2.tc1
            tcenterend = self.before2.tc1 + lenJ1 + lenI2

            srightstart = self.before1.sc2 + lenS + lenJ1 + lenI2 + lenS
            srightend = srightstart + lenJ1 + lenJ2
            trightstart = self.before1.tc2 + lenS + lenJ1 + lenI2 + lenS
            trightend = trightstart + lenJ1 + lenJ2

        if self.direction == LEFT:
            scenterstart = self.before2.sc1
            scenterend = self.before2.sc1 + lenJ1 + lenJ1
            tcenterstart = self.before2.tc1
            tcenterend = self.before2.tc1 + lenJ1 + lenJ1

            srightstart = self.before1.sc2 + lenS + lenJ1 + lenJ1 + lenS
            srightend = srightstart + lenI2 + lenJ2
            trightstart = self.before1.tc2 + lenS + lenJ1 + lenJ1 + lenS
            trightend = trightstart + lenI2 + lenJ2
    
        position = self.before1.position
        self.after2 = Interval(tcenterstart, tcenterend,
                               scenterstart, scenterend, position+1, 'I')
        self.after3 = Interval(trightstart, trightend,
                               srightstart, srightend, position+2, 'I')
