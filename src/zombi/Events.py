from .Interval import Interval

# Types:
T_EVENT = str
T_DIR = bool

# Directions:
RIGHT = False
LEFT = True

# Event types:
TDUP = 'D'  #: Tandem Duplication
FER = 'T'   #: Transfer
LOSS = 'L'  #: Loss
INV = 'I'   #: Inversion
POS = 'P'   #: Transposition
ORIG = 'O'  #: Origination

class GenomeEvent:
    """
    An rearrangement event. Meant to be used as a base class.

    Attributes
    ----------
    etype: str
        the type of event from {TDUP, FER, LOSS, INV, POS, ORIG}
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

    NOTES
    -----
        `before1` and 'sc1' are always considered to be to the left of `before2`
        and `sc2`, so if the event wraps around the indices of `before1` will
        be greater than the indices of `before2`.
    """
    def __init__(self, int1: Interval, int2: Interval, sc1: int, sc2: int,
                 direction: T_DIR, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.before1: Interval = int1
        self.before2: Interval = int2
        self.sc1: int = sc1
        self.sc2: int = sc2
        self.direction: T_DIR = direction

    def wrapsLeft(self) -> bool:
        """
        Does this event wrap left (the `direction` is LEFT and `before1` occurs
        to the left of `before2`)?

        NOTES
        -----
            Assumes no intergenes that wrap from the end to beginning.
        """
        return self.direction == LEFT and self.before1.tc1 < self.before2.tc1

    def wrapsRight(self) -> bool:
        """
        Does this event wrap right (the `direction` is RIGHT and `before1`
        occurs to the right of `before2`)?

        NOTES
        -----
            Assumes no intergenes that wrap from the end to beginning.
        """
        return self.direction == RIGHT and self.before1.tc1 > self.before2.tc1

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

        I is split at `c1` into I0 I1 and J is split at `c2` into J0 J1. Then
        we get

            I0 I1 S J0 J1

        and the inversion produces

            I0 -J0 -S -I1 J1
        """
        lenI0 = self.sc1 - self.before1.sc1
        lenI1 = self.before1.sc2 - self.sc1
        lenJ0 = self.sc2 - self.before2.sc1
        lenJ1 = self.before2.sc2 - self.sc2

        sleftstart = self.before1.sc1
        sleftend = sleftstart + lenI0 + lenJ0
        tleftstart = self.before1.tc1
        tleftend = tleftstart + lenI0 + lenJ0

        srightstart = self.before2.sc1
        srightend = srightstart + lenI1 + lenJ1
        trightstart = self.before2.tc1
        trightend = trightstart + lenI1 + lenJ1

        position = self.before1.position
        self.after1 = Interval(tleftstart, tleftend,
                               sleftstart, sleftend, position, 'I')
        self.after2 = Interval(trightstart, trightend,
                               srightstart, srightend, position+1, 'I')

class TandemDup(EventTwoCuts):
    """
    A tandem duplication event. See the description of `setAfter()` for details.
    
    Attributes
    ----------
    """
    def __init__(self, int1: Interval, int2: Interval, sc1: int, sc2: int,
                 direction: T_DIR, lineage: str, time: float):
        """
        Create a TandemDuplication event.  When instantiating this you must
        provide the arguments for EventTwoCuts and GenomeEvent.

        Parameters
        ----------
        int1: Interval
            the first intergenic interval to be cut
        int2: Interval
            the second intergenic interval to be cut
        sc1: int
            the first cut
        sc2: int
            the second cut
        direction: T_DIR
            the direction, if RIGHT then the second interval is to the right of
            the first
        lineage: str
            the lineage on which the event happened (pendant node name)
        time: float
            the time at which it happened
        """
        super().__init__(int1, int2, sc1, sc2, direction, TDUP, lineage, time)
        self.after1: Interval = None
        self.after2: Interval = None
        self.after3: Interval = None

        self.setAfter()


    def setAfter(self):
        """
        Set the three intergenic regions that exist after the tandem
        duplication.
        If `direction` is RIGHT then consider intergenic regions I = `before1`
        and J = `before2`.  If `direction` is LEFT then intergenic region I
        is `before2` and J is `before1`. I and J are on either side of segment
        S composed of genes and intergenes:

            I S J

        I is split at `c1` into I0 I1 and J is split at `c2` into J0 J1. Then
        we get

            I0 I1 S J0 J1

        and the tandem duplication produces

            I0 I1 S J0 I1 S J0 J1
        """
        if self.direction == RIGHT:
            I = self.before1
            J = self.before2
        else:
            I = self.before2
            J = self.before1

        self.after1 = self.before1

        lenI1 = self.before1.sc2 - self.sc1
        lenS = self.before2.sc1 - self.before1.sc2
        lenJ0 = self.sc2 - self.before2.sc1
        lenJ1 = self.before2.sc2 - self.sc2
        if self.direction == RIGHT:
            scenterstart = self.before2.sc1
            scenterend = self.before2.sc1 + lenJ0 + lenI1
            tcenterstart = self.before2.tc1
            tcenterend = self.before2.tc1 + lenJ0 + lenI1

            srightstart = self.before1.sc2 + lenS + lenJ0 + lenI1 + lenS
            srightend = srightstart + lenJ0 + lenJ1
            trightstart = self.before1.tc2 + lenS + lenJ0 + lenI1 + lenS
            trightend = trightstart + lenJ0 + lenJ1

        if self.direction == LEFT:
            scenterstart = self.before2.sc1
            scenterend = self.before2.sc1 + lenJ0 + lenJ0
            tcenterstart = self.before2.tc1
            tcenterend = self.before2.tc1 + lenJ0 + lenJ0

            srightstart = self.before1.sc2 + lenS + lenJ0 + lenJ0 + lenS
            srightend = srightstart + lenI1 + lenJ1
            trightstart = self.before1.tc2 + lenS + lenJ0 + lenJ0 + lenS
            trightend = trightstart + lenI1 + lenJ1
    
        position = self.before1.position
        self.after2 = Interval(tcenterstart, tcenterend,
                               scenterstart, scenterend, position+1, 'I')
        self.after3 = Interval(trightstart, trightend,
                               srightstart, srightend, position+2, 'I')
        
    def getTotalStr(self):
        return f'{self.before1.tc1, self.before1.tc2} S ' + \
               f'{self.before2.tc1, self.before2.tc2} ->' + \
               f'{self.after1.tc1, self.after1.tc2} S ' + \
               f'{self.after2.tc1, self.after2.tc2} S ' + \
               f'{self.after3.tc1, self.after3.tc2}'

    def getSpecificStr(self):
        return f'{self.before1.sc1, self.before1.sc2} S ' + \
               f'{self.before2.sc1, self.before2.sc2} ->' + \
               f'{self.after1.sc1, self.after1.sc2} S ' + \
               f'{self.after2.sc1, self.after2.sc2} S ' + \
               f'{self.after3.sc1, self.after3.sc2}'