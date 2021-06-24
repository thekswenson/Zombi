import copy
import abc

from Bio.SeqFeature import BeforePosition
from .Interval import Interval

# Types:
T_EVENT = str

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

    @abc.abstractmethod
    def afterToBeforeS(self, sc: int) -> int:
        """
        Given a specific breakpoint coordinate after this event, return the
        same one before.
        """
        pass

    @abc.abstractmethod
    def afterToBeforeT(self, sc: int) -> int:
        """
        Given a specific breakpoint coordinate after this event, return the
        same one before.
        """
        pass

class EventTwoCuts(GenomeEvent):
    """
    An event with two cuts. Meant to be used as a base class.

    NOTES
    -----
        `before1` and 'sc1' are always considered to be to the left of `before2`
        and `sc2`, so if the event wraps around the indices of `before1` will
        be greater than the indices of `before2`.

    ATTRIBUTES
    ----------
    beforeL: Interval
        the first intergenic interval to be cut (Leftmost unless wraps)
    beforeR: Interval
        the second intergenic interval to be cut (Rightmost unless wraps)
    scL: int
        the first cut
    scR: int
        the second cut
    twraplen: int
        the total length of the chromosome
    swraplen: int
        the intergene specific length of the chromosome
    """
    def __init__(self, int1: Interval, int2: Interval, sc1: int, sc2: int,
                 swraplen: int, twraplen: int, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.beforeL: Interval = int1
        self.beforeR: Interval = int2
        self.scL: int = sc1
        self.scR: int = sc2
        self.swraplen: int = swraplen   #: specific wrap length
        self.twraplen: int = twraplen   #: total wrap length

    def wraps(self) -> bool:
        """
        Does this event wrap right (`before1` occurs after `before2`)?

        NOTES
        -----
            Assumes no intergenes that wrap from the end to beginning.
        """
        return self.beforeL.tc1 > self.beforeR.tc1

class Origination(GenomeEvent):
    pass

class Inversion(EventTwoCuts):
    """
    An Inversion event.
    """
    def __init__(self, int1: Interval, int2: Interval, sc1: int, sc2: int,
                 swraplen: int, twraplen: int, sfirstlen:int, tfirstlen:int,
                 ssecondlen: int, tsecondlen: int, lineage: str, time: float):
        """
        Create an Inversion event.  When instantiating this you must also
        provide the arguments for EventTwoCuts and GenomeEvent.

        Parameters
        ----------
        int1 : Interval
            the intergene to be broken at the left
        int2 : Interval
            the intergene to be broken at the right
        sc1 : int
            the specific coordinate where the first intergene is to be cut
        sc2 : int
            the specific coordinate where the second intergene is to be cut
        twraplen: int
            the total length of the chromosome
        swraplen: int
            the intergene specific length of the chromosome
        sfirstlen: int
            if this inversion wraps, then this is the specific length of the
            gene and intergenes that will end up at the end of the genome
            after the inversion.
            (retrieved using `CircularChromosome.inversion_wrap_lengths()`)
        tfirstlen: int
            if this inversion wraps, then this is the total length of the
            gene and intergenes that will end up at the end of the genome
            after the inversion.
        ssecondlen: int
            if this inversion wraps, then this is the specific length of the
            gene and intergenes that will end up at the begining of the genome
            after the inversion.
        tsecondlen: int
            if this inversion wraps, then this is the total length of the
            gene and intergenes that will end up at the beginning of the genome
            after the inversion.
        lineage: str
            the lineage on which the event happened (pendant node name)
        time: float
            the time at which it happened
        """
        super().__init__(int1, int2, sc1, sc2, swraplen, twraplen, INV,
                         lineage, time)
        self.afterL: Interval = None
        self.afterR: Interval = None
        self.sfirstlen = sfirstlen
        self.tfirstlen = tfirstlen
        self.ssecondlen = ssecondlen
        self.tsecondlen = tsecondlen

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
        lenI0 = self.scL - self.beforeL.sc1
        lenI1 = self.beforeL.sc2 - self.scL
        self.lenI1 = lenI1
        lenJ0 = self.scR - self.beforeR.sc1
        self.lenJ0 = lenJ0
        lenJ1 = self.beforeR.sc2 - self.scR
            #lenSs will the be number of intergene nucleotides plus the number
            #of genes in the region S:
        if self.wraps():
            lenSs = (self.swraplen - self.beforeL.sc2) + self.beforeR.sc1 + 1
            lenSt = (self.twraplen - self.beforeL.tc2) + self.beforeR.tc1
        else:
            lenSs = self.beforeR.sc1 - self.beforeL.sc2
            lenSt = self.beforeR.tc1 - self.beforeL.tc2

        if self.wraps():    # [secondlen] -I1 J1 ... I0 -J0 [firstlen]
            self.sshift = self.ssecondlen + lenI1 - self.beforeR.s_breakpoint
            self.tshift = self.tsecondlen + lenI1 - self.beforeR.t_breakpoint

            sleftstart = self.beforeL.sc1 + self.sshift
            sleftend = sleftstart + lenI0 + lenJ0
            tleftstart = self.beforeL.tc1 + self.tshift
            tleftend = tleftstart + lenI0 + lenJ0

            srightstart = self.ssecondlen
            srightend = srightstart + lenI1 + lenJ1
            trightstart = self.tsecondlen
            trightend = trightstart + lenI1 + lenJ1

            sbreakL = self.beforeL.s_breakpoint + self.sshift
            tbreakL = self.beforeL.t_breakpoint + self.tshift
            sbreakR = self.beforeR.s_breakpoint + self.sshift
            tbreakR = self.beforeR.t_breakpoint + self.tshift
        else:
            sleftstart = self.beforeL.sc1
            sleftend = sleftstart + lenI0 + lenJ0
            tleftstart = self.beforeL.tc1
            tleftend = tleftstart + lenI0 + lenJ0

            lenSs = self.beforeR.sc1 - self.beforeL.sc2
            lenSt = self.beforeR.tc1 - self.beforeL.tc2

            srightstart = sleftend + lenSs
            srightend = srightstart + lenI1 + lenJ1
            trightstart = tleftend + lenSt
            trightend = trightstart + lenI1 + lenJ1

            sbreakL = self.beforeL.s_breakpoint
            tbreakL = self.beforeL.t_breakpoint
            sbreakR = self.beforeR.s_breakpoint
            tbreakR = self.beforeR.t_breakpoint

        self.afterL = Interval(tleftstart, tleftend, sleftstart, sleftend,
                               self.beforeL.position, INV, tbreakL, sbreakL)
        self.afterR = Interval(trightstart, trightend, srightstart, srightend,
                               self.beforeR.position, INV, tbreakR, sbreakR)

    def afterToBeforeS(self, sc: int) -> int:
        """
        Given a specific breakpoint coordinate after this Inversion, return the
        same one before.

            I0 I1 S J0 J1 became
            I0 -J0 -S -I1 J1

            or

            S1 J0 J1 ... I0 I1 S0 became
            -S2 -I1 J1 ... I0 -J0 -S3  (e.g. -S2 could have parts of S0 and S1)
        """
        if self.wraps():
            if sc > self.afterL.s_breakpoint:   # sc in -J0 -S3
                right_of_bp = sc - self.afterL.s_breakpoint
                if right_of_bp <= self.beforeR.s_breakpoint:
                    return self.beforeR.s_breakpoint - right_of_bp
                else:
                    return self.swraplen - (right_of_bp - self.beforeR.s_breakpoint - 1)

            elif sc >= self.afterR.s_breakpoint:
                return sc -self.sshift

            else:                               # sc in -S2 -I1
                left_of_bp = self.afterR.s_breakpoint - sc
                if left_of_bp <= self.swraplen - self.beforeL.s_breakpoint:
                    return self.beforeL.s_breakpoint + left_of_bp
                else:
                    return left_of_bp - (self.swraplen - self.beforeL.s_breakpoint) - 1

        else:
            if sc <= self.beforeL.s_breakpoint or sc >= self.beforeR.s_breakpoint:
                return sc
            else:
                return self.beforeL.s_breakpoint + (self.beforeR.s_breakpoint - sc)

    def afterToBeforeT(self, tc: int) -> int:
        """
        Given a total breakpoint coordinate after this Inversion, return the
        same one before.

            I0 I1 S J0 J1 became
            I0 -J0 -S -I1 J1

            or

            S1 J0 J1 ... I0 I1 S0 became
            -S2 -I1 J1 ... I0 -J0 -S3  (e.g. -S2 could have parts of S0 and S1)
        """
        if self.wraps():
            if tc > self.afterL.t_breakpoint:   # tc in -J0 -S3
                right_of_bp = tc - self.afterL.t_breakpoint
                if right_of_bp <= self.beforeR.t_breakpoint:
                    return self.beforeR.t_breakpoint - right_of_bp
                else:
                    return self.twraplen - (right_of_bp - self.beforeR.t_breakpoint)

            elif tc >= self.afterR.t_breakpoint:
                return tc -self.tshift

            else:                               # tc in -S2 -I1
                left_of_bp = self.afterR.t_breakpoint - tc
                if left_of_bp <= self.twraplen - self.beforeL.t_breakpoint:
                    return self.beforeL.t_breakpoint + left_of_bp
                else:
                    return left_of_bp - (self.twraplen - self.beforeL.t_breakpoint)

        else:
            if tc <= self.beforeL.t_breakpoint or tc >= self.beforeR.t_breakpoint:
                return tc
            else:
                return self.beforeL.t_breakpoint + (self.beforeR.t_breakpoint - tc)

class TandemDup(EventTwoCuts):
    """
    A tandem duplication event. See the description of `setAfter()` for details.
    
    ATTRIBUTES
    ----------
    afterL: Interval
        the first intergenic interval after the event (I0 I1, see `setAfter()`)
    afterC: Interval
        the second intergenic interval after the event (J0 I1)
    afterR: Interval
        the third intergenic interval after the event (J0 J1)
    """
    def __init__(self, int1: Interval, int2: Interval, sc1: int, sc2: int,
                 swraplen: int, twraplen: int, lineage: str, time: float):
        """
        Create a TandemDuplication event.  When instantiating this you must
        provide the arguments for EventTwoCuts and GenomeEvent.

        Parameters
        ----------
        int1: Interval
            the left intergenic interval to be cut
        int2: Interval
            the second intergenic interval to be cut
        sc1: int
            the first cut
        sc2: int
            the second cut
        twraplen: int
            the total length of the chromosome
        swraplen: int
            the intergene specific length of the chromosome
        lineage: str
            the lineage on which the event happened (pendant node name)
        time: float
            the time at which it happened
        """
        super().__init__(int1, int2, sc1, sc2, swraplen, twraplen, TDUP,
                         lineage, time)
        self.afterL: Interval = None
        self.afterC: Interval = None
        self.afterR: Interval = None

        self.setAfter()


    def setAfter(self):
        """
        Set the three intergenic regions that exist after the tandem
        duplication.
        Consider intergenic regions I = `before1` and J = `before2`. I and J
        are on either side of segment S composed of genes and intergenes:

            I S J

        I is split at `c1` into I0 I1 and J is split at `c2` into J0 J1. Then
        we get

            I0 I1 S J0 J1

        and the tandem duplication produces

            I0 I1 S J0 I1 S J0 J1
        """
        self.afterL = copy.deepcopy(self.beforeL)

        lenI1 = self.beforeL.sc2 - self.scL #number of (intergene) muclotides
        lenJ0 = self.scR - self.beforeR.sc1 #number of (intergene) muclotides
        lenJ1 = self.beforeR.sc2 - self.scR #number of (intergene) muclotides
            #lenSs will the be number of intergene nucleotides plus the number
            #of genes in the region S:
        if self.wraps():    #beforeL is not to the left if we wrap
            self.lenSs = (self.swraplen - self.beforeL.sc2) + self.beforeR.sc1 + 1
            self.lenSt = (self.twraplen - self.beforeL.tc2) + self.beforeR.tc1
        else:
            self.lenSs = self.beforeR.sc1 - self.beforeL.sc2
            self.lenSt = self.beforeR.tc1 - self.beforeL.tc2

        if self.wraps():
            self.afterL.sc1 += lenI1 + self.lenSs + lenJ0
            self.afterL.sc2 += lenI1 + self.lenSs + lenJ0
            self.afterL.tc1 += lenI1 + self.lenSt + lenJ0
            self.afterL.tc2 += lenI1 + self.lenSt + lenJ0

        scenterstart = self.beforeR.sc1
        scenterend = scenterstart + lenJ0 + lenI1
        scenterbreak = scenterstart + lenJ0
        tcenterstart = self.beforeR.tc1
        tcenterend = tcenterstart + lenJ0 + lenI1
        tcenterbreak = tcenterstart + lenJ0

        srightstart = scenterend + self.lenSs
        srightend = srightstart + lenJ0 + lenJ1
        trightstart = tcenterend + self.lenSt
        trightend = trightstart + lenJ0 + lenJ1
    
        position = self.beforeL.position
        self.afterC = Interval(tcenterstart, tcenterend,
                               scenterstart, scenterend, position+1, 'I',
                               tcenterbreak, scenterbreak)
        self.afterR = Interval(trightstart, trightend,
                               srightstart, srightend, position+2, 'I')

    def afterToBeforeS(self, sc: int) -> int:
        """
        Given a specific breakpoint coordinate after this TandemDup, return the
        same one before.

            I0 I1 S J0 J1 became
            I0 I1 S J0 I1 S J0 J1

        The breakpoint between J0 and I1 is the only ambiguous breakpoint. We
        arbitarily map to the left breakpoint I0 I1 (rather than J0 J1).
        """
        if self.afterL.inSpecific(sc):                  # I0 I1
            return self.beforeL.sc1 + (sc - self.afterL.sc1)
        elif self.afterC.inSpecific(sc):
            if sc < self.afterC.s_breakpoint:           # J0
                return self.beforeR.sc1 + (sc - self.afterC.sc1)
            else:                                       # I1
                return self.beforeL.sc2 - (self.afterC.sc2 - sc)
        elif self.afterR.inSpecific(sc):                # J0 J1
            return self.beforeR.sc1 + (sc - self.afterR.sc1)
        elif sc > self.afterC.sc2:                      # to the right
            return sc - self.lenSs - self.afterC.specificLen()
        else:                                           # to the left
            return sc

    def afterToBeforeT(self, tc: int) -> int:
        """
        Given a total breakpoint coordinate after this TandemDup, return the
        same one before.

            I0 I1 S J0 J1 became
            I0 I1 S J0 I1 S J0 J1

        The breakpoint between J0 and I1 is the only ambiguous breakpoint. We
        arbitarily map to the left breakpoint I0 I1 (rather than J0 J1).
        """
        if self.afterL.inTotal(tc):             # I0 I1
            return self.beforeL.tc1 + (tc - self.afterL.tc1)
        elif self.afterC.inTotal(tc):
            if tc < self.afterC.t_breakpoint:   # J0
                return self.beforeR.tc1 + (tc - self.afterC.tc1)
            else:                               # I1
                return self.beforeL.tc2 - (self.afterC.tc2 - tc)
        elif self.afterR.inTotal(tc):           # J0 J1
            return self.beforeR.tc1 + (tc - self.afterR.tc1)
        elif tc > self.afterC.tc2:              # to the right
            return tc - self.lenSt - self.afterC.totalLen()
        else:                                   # to the left
            return tc

    def getTotalStr(self):
        return f'{self.beforeL.tc1, self.beforeL.tc2} S ' + \
               f'{self.beforeR.tc1, self.beforeR.tc2} -> ' + \
               f'{self.afterL.tc1, self.afterL.tc2} S ' + \
               f'{self.afterC.tc1, self.afterC.tc2} S ' + \
               f'{self.afterR.tc1, self.afterR.tc2}'

    def getSpecificStr(self):
        return f'{self.beforeL.sc1, self.beforeL.sc2} S ' + \
               f'{self.beforeR.sc1, self.beforeR.sc2} -> ' + \
               f'{self.afterL.sc1, self.afterL.sc2} S ' + \
               f'{self.afterC.sc1, self.afterC.sc2} S ' + \
               f'{self.afterR.sc1, self.afterR.sc2}'

    def __repr__(self):
        return f'{repr(self.beforeL)} {repr(self.beforeR)} {self.twraplen} ' +\
               f'{self.swraplen} {self.twraplen} {self.lineage} {self.time}'
