import copy
import abc
from os import initgroups
from typing import List

from zombi.Genomes import Intergene

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

class EventOneCut(GenomeEvent):
    """
    An event with one cut. Meant to be used as a base class.

    ATTRIBUTES
    ----------
    before: Interval
        the intergenic interval to be cut
    """
    def __init__(self, interval: Interval, sbp: int, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert interval.inSpecific(sbp)
        self.before: Interval = interval
        self.sbp: int = sbp

class EventTwoCuts(GenomeEvent):
    """
    An event with two cuts. Meant to be used as a base class.

    NOTES
    -----
        `before1` and 'sbp1' are always considered to be to the left of
        `before2` and `sbp2`, so if the event wraps around the indices of
        `before1` will be greater than the indices of `before2`.

    ATTRIBUTES
    ----------
    beforeL: Interval
        the first intergenic interval to be cut (Leftmost unless wraps). This
        contains information about the interval and breakpoint within it.
    beforeR: Interval
        the second intergenic interval to be cut (Rightmost unless wraps)
    sbpL: int
        the first (specific) cut
    sbpR: int
        the second (specific) cut
    tbpL: int
        the first (total) cut
    tbpR: int
        the second (total) cut
    twraplen: int
        the total length of the chromosome
    swraplen: int
        the intergene specific length of the chromosome
    """
    def __init__(self, int1: Interval, int2: Interval, sbp1: int, sbp2: int,
                 swraplen: int, twraplen: int, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.beforeL: Interval = int1
        self.beforeR: Interval = int2
        self.sbpL: int = sbp1           #: specific break-point on Left
        self.sbpR: int = sbp2           #: specific break-point on Right
        self.tbpL: int = int1.tc1 + (sbp1 - int1.sc1)
        self.tbpR: int = int2.tc1 + (sbp2 - int2.sc1)
        self.swraplen: int = swraplen   #: specific wrap length
        self.twraplen: int = twraplen   #: total wrap length

    def wraps(self) -> bool:
        """
        Does this event wrap around to the right (`before1` occurs after
        `before2`)?

        NOTES
        -----
            Assumes no intergenes that wrap from the end to beginning.
        """
        return self.beforeL.tc1 > self.beforeR.tc1

    def assertS(self, sc: int) -> int:
        """
        Return the given value after checking it is within the bounds.
        """
        assert 0 <= sc <= self.swraplen
        return sc

    def assertT(self, tc: int) -> int:
        """
        Return the given value after checking it is within the bounds.
        """
        assert 0 <= tc <= self.twraplen
        return tc

#-- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - --
class Origination(EventOneCut):
    """
    An origination, where a new gene in placed within an intergene.

    Attributes
    ----------
    afterL: Interval
        the left intergenic interval after the cut
    afterR: Interval
        the right intergenic interval after the cut
    """
    def __init__(self, interval: Interval, sbp: int, genelen:int,
                 lineage: str, time: float):
        """
        Create an Origination event.

        Parameters
        ----------
        interval: Interval
            the intergenic interval to put the new gene in
        sbp: int
            the exact spot to put it cut
        genelen: int
            the length of the inserted gene
        lineage: str
            the lineage on which the event happened (pendant node name)
        time: float
            the time at which it happened
        """
        super().__init__(interval, sbp, ORIG, lineage, time)
        self.afterL: Interval = None
        self.afterR: Interval = None
        self.genelen = genelen

        self.setAfter()

    def setAfter(self):
        """
        Compute the resulting intergenic intervals after the cut.

        Consider intergenic regions I = `before`:

            ... I ...

        I is split at `self.sbp` into I0 I1.  And the origination produces:

            ... I0 G I1 ...
        """
        lenI0 = self.sbp - self.before.sc1
        lenI1 = self.before.sc2 - self.sbp
        self.afterL = Interval(self.before.tc1, self.before.tc1 + lenI0,
                               self.before.sc1, self.sbp,
                               self.before.position, 'I')

        tstart = (self.before.t_bp - lenI1) + self.genelen + 1
        self.afterR = Interval(tstart, tstart + lenI1,
                               self.sbp+1, self.sbp+1 + lenI1,
                               self.before.position+1, 'I')

    def afterToBeforeS(self, sc: int) -> int:
        """
        Given a specific breakpoint coordinate after this Origination, return
        the same one before.

            ... I ...       became
            ... I0 G I1 ...
        """
        if sc > self.afterL.sc2:
            return sc-1
        
        return sc

    def afterToBeforeT(self, tc: int) -> int:
        """
        Given a total breakpoint coordinate after this Origination, return
        the same one before.

            ... I ...       became
            ... I0 G I1 ...
        """
        if self.afterL.tc2 < tc < self.afterR.tc1:
            raise(MapOriginError('tried to map a gene coordinate'))
        elif tc >= self.afterR.tc1:
            return tc - self.genelen
        
        return tc

#-- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - --
class Loss(EventTwoCuts):
    """
    A loss event.

    Attributes
    ----------
    after: Interval
        the intergenic interval after the cut
    pseudogenize: bool
        the region between the breakpoints was pseudogenized by this Loss
    pseudolist: List[interval]
        list of intergenes in the region to be pseudogenized
    genebps: int
        number of gene breakpoints in the pseudogenized region
    after_sbptL: int
        the specific breakpoint (on the left) after a pseudogenization
    after_sbptR: int
        the specific breakpoint (on the right) after a pseudogenization
    after_tbptL: int
        the total breakpoint (on the left) after a pseudogenization
    after_tbptR: int
        the total breakpoint (on the right) after a pseudogenization
    """
    def __init__(self, int1: Interval, int2: Interval, sbp1: int, sbp2: int,
                 swraplen: int, twraplen: int, lineage: str, time: float,
                 pseudogenize=False, pseudolist: List[Intergene]=None):
        """
        Create a Loss event. Either cut out everything between `sbp1` and
        `sbp2`, or turn everything in that region into a big intergene,
        depending on the value `pseudogenize`.

        Parameters
        ----------
        int1: Interval
            the left intergenic interval to be cut
        int2: Interval
            the second intergenic interval to be cut
        sbp1: int
            the first cut
        sbp2: int
            the second cut
        pseudogenize: bool
            turn everything into one big intergene
        pseudolist: List[Interval]
            the intergenes to pseudogenize
        twraplen: int
            the total length of the chromosome
        swraplen: int
            the intergene specific length of the chromosome
        lineage: str
            the lineage on which the event happened (pendant node name)
        time: float
            the time at which it happened
        """
        super().__init__(int1, int2, sbp1, sbp2, swraplen, twraplen, LOSS,
                         lineage, time)
        self.after: Interval = None

        self.pseudogenize = pseudogenize
        self.pseudolist = copy.deepcopy(pseudolist)
        assert not pseudolist or pseudogenize, 'psuedolist given but pseudogenize is False'
        self.genebps: int = -1  #: num gene breakpoints in pseudogenized region

        self.after_sbpL: int = -1
        self.after_sbpR: int = -1
        self.after_tbpL: int = -1
        self.after_tbpR: int = -1

        self.setAfter()

    def setAfter(self):
        """
        Compute the resulting intergenic interval after the cut.

        Consider intergenic regions I = `before1` and J = `before2` on either
        side of segment S:

            I S J

        I is split at `sbpL` into I0 I1 and J is split at `sbpR` into J0 J1.
        Then we have

            I0 I1 S J0 J1

        and the loss produces

            I0 J1
        """
        lenI0 = self.sbpL - self.beforeL.sc1
        lenJ1 = self.beforeR.sc2 - self.sbpR

        sstart = self.beforeL.sc1
        tstart = self.beforeL.tc1

        position = self.beforeL.position
        newsbp = self.sbpL
        newtbp = self.tbpL
        if self.pseudogenize:
            if self.wraps():
                position -= self.beforeR.position + 1

                sstart -= self.beforeR.sc2 + 1
                tstart -= self.beforeR.tc2
                send = sstart + self.beforeR.tc2 + (self.twraplen - self.beforeL.tc1)
                tend = self.twraplen

                newsbp -= self.beforeR.sc2 + 1
                newtbp -= self.beforeR.tc2

                self.after_sbpL = self.sbpL - self.beforeR.sc2 - 1
                self.after_sbpR = send - (self.beforeR.sc2 - self.sbpR)
                self.after_tbpL = self.tbpL - self.beforeR.tc2
                self.after_tbpR = tend - (self.beforeR.tc2 - self.tbpR)
            else:
                send = sstart + (self.beforeR.tc2 - tstart)
                tend = tstart + (self.beforeR.tc2 - tstart)

                self.genebps = ((self.beforeR.tc1 - self.beforeL.tc2) -
                                (self.beforeR.sc1 - self.beforeL.sc2))

                self.after_sbpL = self.sbpL
                self.after_sbpR = self.sbpR + self.genebps
                self.after_tbpL = self.tbpL
                self.after_tbpR = self.tbpR
        else:
            if self.wraps():
                position -= self.beforeR.position + 1
                newsbp -= self.beforeR.sc2 + 1
                newtbp -= self.beforeR.tc2
                sstart -= self.beforeR.sc2 + 1
                tstart -= self.beforeR.tc2

            send = sstart + lenI0 + lenJ1
            tend = tstart + lenI0 + lenJ1

        self.after = Interval(tstart, tend, sstart, send,
                              position, 'I', newtbp, newsbp)

    def afterToBeforeS(self, sc: int) -> int:
        """
        Given a specific breakpoint coordinate after this Loss, return the
        same one before.

            .. I0 I1 S J0 J1 .. became
            .. I0 J1 ..

            or

            S1 J0 J1 .. I0 I1 S0 became
            .. I0 J1

        (with no pseudogenization, otherwise it's different)
        """
        if self.pseudogenize:
            if self.wraps():
                if sc <= self.after_sbpL:                   # in I0 or before
                    return self.assertS(sc + self.beforeR.sc2 + 1)
                
                elif sc >= self.after_sbpR:                 # in J1
                    return self.assertS(self.beforeR.s_bp + (sc - self.after_sbpR))
                elif (retval := self.inPseudoIntergene(sc)) >= 0:
                    return self.assertS(retval)             # in I1 SO or S1 J0
                else:
                    raise(MapPseudogeneError(f'specific coordinate {sc} inside pseudogenized region'))
            else:
                if sc >= self.after_sbpR:                   # in J1 or to right
                    return self.assertS(sc - self.genebps)

                elif sc <= self.after_sbpL:                 # in I0 or to left
                    return self.assertS(sc)
                elif (retval := self.inPseudoIntergene(sc)) >= 0:
                    return self.assertS(retval)             # in I1 S J0
                else:
                    raise(MapPseudogeneError(f'specific coordinate {sc} inside pseudogenized region'))
        else:
            if self.wraps():
                if sc > self.after.s_bp:                # in J1
                    return self.assertS(self.beforeR.s_bp + (sc - self.after.s_bp))
                else:                                   # in I0 or to left
                    return self.assertS(sc + self.beforeR.sc2 + 1)
            else:
                if sc > self.after.s_bp:                # in J1 or to right
                    return self.assertS(sc + self.beforeR.s_bp - self.beforeL.s_bp)
                else:                                   # in I0 or to left
                    return self.assertS(sc)

    def inPseudoIntergene(self, sc: int) -> int:
        """
        If the given specific breakpoint coordinate maps to an intergene that
        was psuedogenizes, then return the coordinate before the loss.

        Parameters
        ----------
        sc : int
            the specific coordinate to query

        Returns
        -------
        int
            the specific coordinate before the Loss
        """
        assert self.pseudogenize, 'Loss was not a pseudogenization'

        if self.wraps():
            assert(self.after_sbpL <= sc <= self.after_sbpR)
                # in I1
            if sc <= self.after_sbpL + (self.beforeL.sc2 - self.beforeL.s_bp):
                return self.assertS(self.beforeL.s_bp + (sc - self.after_sbpL))

                # in J0
            if sc >= self.after_sbpR - (self.beforeR.s_bp - self.beforeR.sc1):
                return self.assertS(self.beforeR.s_bp - (self.after_sbpR - sc))

                # in an intergene of S0 or S1
            tcbefore = self.beforeL.t_bp + sc - self.after_sbpL
            if tcbefore > self.twraplen:
                tcbefore -= self.twraplen

            for intergene in self.pseudolist:
                if intergene.inTotal(tcbefore):
                    return self.assertS(intergene.sc1 + (tcbefore - intergene.tc1))
        else:
            assert self.after_sbpL <= sc <= self.after_sbpR

            tc = self.after_tbpL + (sc - self.after_sbpL)
            for intergene in self.pseudolist:
                if intergene.inTotal(tc):                   # in S
                    return self.assertS(intergene.sc1 + (tc - intergene.tc1))

            if self.beforeL.inSpecific(sc):                 # in I1
                return self.assertS(sc)
            
            elif self.beforeR.inTotal(tc):
                return self.assertS(sc - self.genebps)      # in J0

        return -1

    def afterToBeforeT(self, tc: int) -> int:
        """
        Given a total breakpoint coordinate after this Loss, return the
        same one before.

            .. I0 I1 S J0 J1 .. became
            .. I0 J1 ..

            or

            S1 J0 J1 .. I0 I1 S0 became
            .. I0 J1

        (with no pseudogenization, otherwise it's different)

        NOTES
        -----
            Any total coordinate mapped within the pseudogenized region will
            raise a NotImplemented error (since the intergenes in I1 S J0 are
            not mapped to anything, as they are in afterToBeforeS())
        """
        if self.pseudogenize:
            if self.wraps():
                if tc <= self.after_tbpL:                   # in I0 or before
                    return self.assertT(tc + self.beforeR.tc2)
                
                elif tc >= self.after_tbpR:                 # in J1
                    return self.assertT(self.beforeR.t_bp + (tc - self.after_tbpR))
                else:
                    raise(NotImplementedError)
                    #raise(MapPseudogeneError(f'total coordinate {tc} inside pseudogenized region'))
            else:
                if tc >= self.after_tbpR or tc <= self.after_tbpL:
                    return self.assertT(tc)
                else:
                    raise(NotImplementedError)
                    #raise(MapPseudogeneError(f'total coordinate {tc} inside pseudogenized region'))
        else:
            if self.wraps():
                if tc > self.after.t_bp:                # in J1
                    return self.assertT(self.beforeR.t_bp + (tc - self.after.t_bp))
                else:                                   # in I0 or to left
                    return self.assertT(tc + self.beforeR.tc2)
            else:
                if tc > self.after.t_bp:                # in J1 or to right
                    return self.assertT(tc + self.beforeR.t_bp - self.beforeL.t_bp)
                else:                                   # in I0 or to left
                    return self.assertT(tc)

class MapPseudogeneError(Exception):
    pass
class MapOriginError(Exception):
    pass

#-- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - --
class Transposition(EventTwoCuts):
    """
    A transposition genome event.

    ATTRIBUTES
    ----------
    beforeH: Interval
        the copied segment will be place here, in this intergenic region
    sbpH: int
        the copied segment will be place here, at this specific breakpoint
    tbpH: int
        the copied segment will be place here, at this total breakpoint
    afterL: Interval
        the intergenic interval to the left of the transposed segment
        (see setAfter())
    afterR: Interval
        the intergenic interval to the right of the transposed segment
    afterH: Interval
        here is where the segment was moved from (I0 J1)
    numpositions: int
        the total number of intergenes in the chromosome
    """
    def __init__(self, int1: Interval, int2: Interval, sbp1: int, sbp2: int,
                 int3: Interval, sbp3: int, numintergenes: int,
                 swraplen: int, twraplen: int, lineage: str, time: float):
        """
        Create a Transposition event.  

        Parameters
        ----------
        int1 : Interval
            the intergene to be broken at the left
        int2 : Interval
            the intergene to be broken at the right
        sbp1 : int
            the specific coordinate where the first intergene is to be cut
        sbp2 : int
            the specific coordinate where the second intergene is to be cut
        int3 : Interval
            the intergene where the segment will be placed
        sbp3 : int
            intergene specific breakpoint coordinate where to place the segment
        numintergenes : int
            the number of intergene positions in the chromosome
        twraplen: int
            the total length of the chromosome
        swraplen: int
            the intergene specific length of the chromosome
        lineage: str
            the lineage on which the event happened (pendant node name)
        time: float
            the time at which it happened
        """
        super().__init__(int1, int2, sbp1, sbp2, swraplen, twraplen, INV,
                         lineage, time)
        self.beforeH: Interval = int3
        self.sbpH = sbp3
        self.tbpH: int = int3.tc1 + (sbp3 - int3.sc1)

        self.afterL: Interval = None
        self.afterR: Interval = None
        self.afterH: Interval = None

        self.numintergenes = numintergenes

        self.setAfter()

    def setAfter(self):
        """
        Set the three intergenic regions that exist after the Transposition.
        Consider intergenic regions I = `before1` and J = `before2` on either
        side of segment S and K = `before3` where the transposed region will
        land.

            I S J ... K

        I is split at `self.sbpL` into I0 I1 and J is split at `self.sbpR` into
        J0 J1. K is split at 
        Then we have

            I0 I1 S J0 J1 ... K0 K1
            (K0 K1 ... I0 I1 S J0 J1)

        and the transposition produces

            I0 J1 ... K0 I1 S J0 K1
            (K0 I1 S J0 K1 ... I0 J1).

        When a transposition wraps around it looks like this

            S1 J0 J1 ... K0 K1 ... I0 I1 S0 becomes
            ... K0 I1 S0 S1 J0 K1 ... I0 J1
        """
        lenI0 = self.sbpL - self.beforeL.sc1
        lenI1 = self.beforeL.sc2 - self.sbpL
        lenJ0 = self.sbpR - self.beforeR.sc1
        lenJ1 = self.beforeR.sc2 - self.sbpR
        lenK0 = self.sbpH - self.beforeH.sc1
        lenK1 = self.beforeH.sc2 - self.sbpH

        Spositions = self.beforeR.position - self.beforeL.position - 1

            #lenSs will the be number of intergene nucleotides plus the number
            #of genes in the region S:
        if self.wraps():
            lenS0s = self.swraplen - self.beforeL.sc2
            lenS1s = self.beforeR.sc1 + 1
            lenSs = lenS0s + lenS1s

            lenS0t = self.twraplen - self.beforeL.tc2
            lenS1t = self.beforeR.tc1
            lenSt =  lenS0t + lenS1t

            self.lenS0s = lenS0s
            self.lenS1s = lenS1s
            self.lenS0t = lenS0t
            self.lenS1t = lenS1t
        else:
            lenSs = self.beforeR.sc1 - self.beforeL.sc2
            lenSt = self.beforeR.tc1 - self.beforeL.tc2

        self.lenI1 = lenI1
        self.lenJ0 = lenJ0
        self.lenJ1 = lenJ1
        self.lenSs = lenSs
        self.lenSt = lenSt

        if self.wraps():
            s_herestart = self.beforeL.sc1 + lenI1 + lenS0s - lenJ1
            s_hereend = s_herestart + lenI0 + lenJ1
            t_herestart = self.beforeL.tc1 + lenI1 + lenS0t - lenJ1
            t_hereend = t_herestart + lenI0 + lenJ1
            assert t_hereend == self.twraplen
            
            s_leftstart = self.beforeH.sc1 - (lenS1s + lenJ0 + lenJ1)
            s_leftend = s_leftstart + (self.beforeH.s_bp - self.beforeH.sc1) + lenI1
            t_leftstart = self.beforeH.tc1 - (lenS1t + lenJ0 + lenJ1)
            t_leftend = t_leftstart + (self.beforeH.t_bp - self.beforeH.tc1) + lenI1

            s_rightstart = self.beforeH.s_bp - lenJ1 + lenI1 + lenS0s - lenJ0
            s_rightend = s_rightstart + (self.beforeH.sc2 - self.beforeH.s_bp) + lenJ0
            t_rightstart = self.beforeH.t_bp - lenJ1 + lenI1 + lenS0t - lenJ0
            t_rightend = t_rightstart + (self.beforeH.tc2 - self.beforeH.t_bp) + lenJ0

            leftposition = 0
            rightposition = self.beforeH.position + (self.numintergenes - self.beforeL.position - 1)
            hereposition = self.numintergenes - 1
        else:
            if self.sbpH > self.sbpR:   # I0 I1 S J0 J1 ... K0 K1
                s_herestart = self.beforeL.sc1
                s_hereend = s_herestart + lenI0 + lenJ1
                t_herestart = self.beforeL.tc1
                t_hereend = t_herestart + lenI0 + lenJ1

                s_leftstart = self.beforeH.sc1 - (lenSs + lenI1 + lenJ0)
                s_leftend = s_leftstart + lenK0 + lenI1
                t_leftstart = self.beforeH.tc1 - (lenSt + lenI1 + lenJ0)
                t_leftend = t_leftstart + lenK0 + lenI1

                leftposition = self.beforeH.position + Spositions
                rightposition = self.beforeR.position
                hereposition = self.beforeL.position

            elif self.sbpH < self.sbpL: # K0 K1 ... I0 I1 S J0 J1
                s_herestart = self.beforeL.sc1 + lenI1 + lenJ0 + lenSs
                s_hereend = s_herestart + lenI0 + lenJ1
                t_herestart = self.beforeL.tc1 + lenI1 + lenJ0 + lenSt
                t_hereend = t_herestart + lenI0 + lenJ1

                s_leftstart = self.beforeH.sc1
                s_leftend = s_leftstart + lenK0 + lenI1
                t_leftstart = self.beforeH.tc1
                t_leftend = t_leftstart + lenK0 + lenI1

                leftposition = self.beforeH.position
                rightposition = leftposition + Spositions
                hereposition = self.beforeR.position

            else:
                raise(Exception('unexpected relationship between bp coordinates'))

            s_rightstart = s_leftend + lenSs
            s_rightend = s_rightstart + lenJ0 + lenK1
            t_rightstart = t_leftend + lenSt
            t_rightend = t_rightstart + lenJ0 + lenK1

        sbreakL = s_leftstart + lenK0
        sbreakR = s_rightstart + lenJ0
        sbreakH = s_herestart + lenI0
        tbreakL = t_leftstart + lenK0
        tbreakR = t_rightstart + lenJ0
        tbreakH = t_herestart + lenI0

        self.afterH = Interval(t_herestart, t_hereend, s_herestart, s_hereend,
                               hereposition, 'I', tbreakH, sbreakH)
        self.afterL = Interval(t_leftstart, t_leftend, s_leftstart, s_leftend,
                               leftposition, 'I', tbreakL, sbreakL)
        self.afterR = Interval(t_rightstart, t_rightend, s_rightstart, s_rightend,
                               rightposition, 'I', tbreakR, sbreakR)

    def afterToBeforeS(self, sc: int) -> int:
        """
        Given a specific breakpoint coordinate after this Transposition, return
        the same one before.

            I0 I1 S J0 J1 ... K0 K1
            I0 J1 ... K0 I1 S J0 K1

            or

            K0 K1 ... I0 I1 S J0 J1
            K0 I1 S J0 K1 ... I0 J1

            or

            S1 J0 J1 ... K0 K1 ... I0 I1 S0
            ... K0 I1 S0 S1 J0 K1 ... I0 J1

        """
        if self.wraps():
            if sc <= self.afterL.s_bp:                              # ... K0
                return sc + self.lenS1s + self.lenJ0 + self.lenJ1
            elif sc <= self.afterL.s_bp + self.lenI1 + self.lenS0s: # I1 S0
                return self.beforeL.s_bp + (sc - self.afterL.s_bp)
            elif sc < self.afterR.s_bp:                             # S1 J0
                return sc - (self.afterL.s_bp + self.lenI1 + self.lenS0s) - 1
            elif sc <= self.afterH.s_bp:                            # K1 ... I0
                return sc - (self.lenI1 + self.lenS0s) + self.lenJ1
            else:                                                   # J1
                return (sc - self.afterH.s_bp) + self.lenS1s + self.lenJ0 - 1
        else:
            if self.sbpH > self.sbpR:               # I0 J1 ... K0 I1 S J0 K1
                if self.afterH.s_bp < sc <= self.afterL.s_bp:       # J1 ... K0
                    return sc + self.lenI1 + self.lenSs + self.lenJ0
                elif self.afterL.s_bp < sc < self.afterR.s_bp:      # I1 S J0
                    return self.afterH.s_bp + (sc - self.afterL.s_bp)
                else:
                    return sc
            else:                                   # K0 I1 S J0 K1 ... I0 J1
                if self.afterL.s_bp < sc < self.afterR.s_bp:        # I1 S J0
                    return self.beforeL.s_bp + (sc - self.afterL.s_bp)
                elif self.afterR.s_bp <= sc < self.afterH.s_bp:     # K1 ... I0
                    return sc - (self.lenI1 + self.lenSs + self.lenJ0)
                else:
                    return sc

    def afterToBeforeT(self, tc: int) -> int:
        """
        Given a total breakpoint coordinate after this Transposition, return the
        same one before.

            I0 I1 S J0 J1 ... K0 K1
            I0 J1 ... K0 I1 S J0 K1

            or

            K0 K1 ... I0 I1 S J0 J1
            K0 I1 S J0 K1 ... I0 J1

            or

            S1 J0 J1 ... K0 K1 ... I0 I1 S0
            ... K0 I1 S0 S1 J0 K1 ... I0 J1
        """
        if self.wraps():
            if tc <= self.afterL.t_bp:                              # ... K0
                return tc + self.lenS1t + self.lenJ0 + self.lenJ1
            elif tc <= self.afterL.t_bp + self.lenI1 + self.lenS0t: # I1 S0
                return self.beforeL.t_bp + (tc - self.afterL.t_bp)
            elif tc < self.afterR.t_bp:                             # S1 J0
                return tc - (self.afterL.t_bp + self.lenI1 + self.lenS0t)
            elif tc <= self.afterH.t_bp:                            # K1 ... I0
                return tc - (self.lenI1 + self.lenS0t) + self.lenJ1
            else:                                                   # J1
                return (tc - self.afterH.t_bp) + self.lenS1t + self.lenJ0
        else:
            if self.tbpH > self.tbpR:               # I0 J1 ... K0 I1 S J0 K1
                if self.afterH.t_bp < tc <= self.afterL.t_bp:       # J1 ... K0
                    return tc + self.lenI1 + self.lenSt + self.lenJ0
                elif self.afterL.t_bp < tc < self.afterR.t_bp:      # I1 S J0
                    return self.afterH.t_bp + (tc - self.afterL.t_bp)
                else:
                    return tc
            else:                                   # K0 I1 S J0 K1 ... I0 J1
                if self.afterL.t_bp < tc < self.afterR.t_bp:        # I1 S J0
                    return self.beforeL.t_bp + (tc - self.afterL.t_bp)
                elif self.afterR.t_bp <= tc < self.afterH.t_bp:     # K1 ... I0
                    return tc - (self.lenI1 + self.lenSt + self.lenJ0)
                else:
                    return tc

#-- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - --
class Inversion(EventTwoCuts):
    """
    An Inversion event.

    ATTRIBUTES
    ----------
    afterL: Interval
        the first intergenic interval after the event (I0 -J0, see `setAfter()`)
    afterR: Interval
        the second intergenic interval after the event (-I1 J1)
    """
    def __init__(self, int1: Interval, int2: Interval, sbp1: int, sbp2: int,
                 swraplen: int, twraplen: int, sfirstlen:int, tfirstlen:int,
                 ssecondlen: int, tsecondlen: int, lineage: str, time: float):
        """
        Create an Inversion event.

        Parameters
        ----------
        int1 : Interval
            the intergene to be broken at the left
        int2 : Interval
            the intergene to be broken at the right
        sbp1 : int
            the specific coordinate where the first intergene is to be cut
        sbp2 : int
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
        super().__init__(int1, int2, sbp1, sbp2, swraplen, twraplen, INV,
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
        side of segment S:

            I S J

        I is split at `sbpL` into I0 I1 and J is split at `sbpR` into J0 J1.
        Then we have

            I0 I1 S J0 J1

        and the inversion produces

            I0 -J0 -S -I1 J1
        """
        lenI0 = self.sbpL - self.beforeL.sc1
        lenI1 = self.beforeL.sc2 - self.sbpL
        lenJ0 = self.sbpR - self.beforeR.sc1
        lenJ1 = self.beforeR.sc2 - self.sbpR
            #lenSs will the be number of intergene nucleotides plus the number
            #of genes in the region S:
        if self.wraps():
            lenSs = (self.swraplen - self.beforeL.sc2) + self.beforeR.sc1 + 1
            lenSt = (self.twraplen - self.beforeL.tc2) + self.beforeR.tc1
        else:
            lenSs = self.beforeR.sc1 - self.beforeL.sc2
            lenSt = self.beforeR.tc1 - self.beforeL.tc2

        if self.wraps():    # [secondlen] -I1 J1 ... I0 -J0 [firstlen]
            self.sshift = self.ssecondlen + lenI1 - self.sbpR
            self.tshift = self.tsecondlen + lenI1 - self.tbpR

            sleftstart = self.beforeL.sc1 + self.sshift
            sleftend = sleftstart + lenI0 + lenJ0
            tleftstart = self.beforeL.tc1 + self.tshift
            tleftend = tleftstart + lenI0 + lenJ0

            srightstart = self.ssecondlen
            srightend = srightstart + lenI1 + lenJ1
            trightstart = self.tsecondlen
            trightend = trightstart + lenI1 + lenJ1

            sbreakL = self.sbpL + self.sshift
            tbreakL = self.tbpL + self.tshift
            sbreakR = self.sbpR + self.sshift
            tbreakR = self.tbpR + self.tshift
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

            sbreakL = self.sbpL
            tbreakL = self.tbpL
            sbreakR = self.sbpR
            tbreakR = self.tbpR

        self.afterL = Interval(tleftstart, tleftend, sleftstart, sleftend,
                               self.beforeL.position, 'I', tbreakL, sbreakL)
        self.afterR = Interval(trightstart, trightend, srightstart, srightend,
                               self.beforeR.position, 'I', tbreakR, sbreakR)

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
            if sc > self.afterL.s_bp:   # sc in -J0 -S3
                right_of_bp = sc - self.afterL.s_bp
                if right_of_bp <= self.sbpR:
                    return self.assertS(self.sbpR - right_of_bp)
                else:
                    return self.assertS(self.swraplen - (right_of_bp - self.sbpR - 1))

            elif sc >= self.afterR.s_bp:
                return self.assertS(sc -self.sshift)

            else:                               # sc in -S2 -I1
                left_of_bp = self.afterR.s_bp - sc
                if left_of_bp <= self.swraplen - self.sbpL:
                    return self.assertS(self.sbpL + left_of_bp)
                else:
                    return self.assertS(left_of_bp - (self.swraplen - self.sbpL) - 1)

        else:
            if sc <= self.sbpL or sc >= self.sbpR:
                return self.assertS(sc)
            else:
                return self.assertS(self.sbpL + (self.sbpR - sc))

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
            if tc > self.afterL.t_bp:   # tc in -J0 -S3
                right_of_bp = tc - self.afterL.t_bp
                if right_of_bp <= self.tbpR:
                    return self.assertT(self.tbpR - right_of_bp)
                else:
                    return self.assertT(self.twraplen - (right_of_bp - self.tbpR))

            elif tc >= self.afterR.t_bp:
                return self.assertT(tc -self.tshift)

            else:                               # tc in -S2 -I1
                left_of_bp = self.afterR.t_bp - tc
                if left_of_bp <= self.twraplen - self.tbpL:
                    return self.assertT(self.tbpL + left_of_bp)
                else:
                    return self.assertT(left_of_bp - (self.twraplen - self.tbpL))

        else:
            if tc <= self.tbpL or tc >= self.tbpR:
                return self.assertT(tc)
            else:
                return self.assertT(self.tbpL + (self.tbpR - tc))

#-- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - --
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
    lenSs: int
        number of intergene-specific bases between the broken intergenes
    lenSt: int
        total number of bases between the broken intergenes
    numintergenes: int
        number of intergenes that are being duplicated
    """
    def __init__(self, int1: Interval, int2: Interval, sbp1: int, sbp2: int,
                 numintergenes: int, swraplen: int, twraplen: int,
                 lineage: str, time: float):
        """
        Create a tandem duplication event.

        Parameters
        ----------
        int1: Interval
            the left intergenic interval to be cut
        int2: Interval
            the second intergenic interval to be cut
        sbp1: int
            the first cut
        sbp2: int
            the second cut
        numintergenes: int
            the number of intergenes that are being duplicated
        twraplen: int
            the total length of the chromosome
        swraplen: int
            the intergene specific length of the chromosome
        lineage: str
            the lineage on which the event happened (pendant node name)
        time: float
            the time at which it happened
        """
        super().__init__(int1, int2, sbp1, sbp2, swraplen, twraplen, TDUP,
                         lineage, time)
        self.afterL: Interval = None
        self.afterC: Interval = None
        self.afterR: Interval = None

        self.numintergenes: int = numintergenes

        self.setAfter()


    def setAfter(self):
        """
        Set the three intergenic regions that exist after the tandem
        duplication.
        Consider intergenic regions I = `before1` and J = `before2`. I and J
        are on either side of segment S composed of genes and intergenes:

            I S J

        I is split at `sbpL` into I0 I1 and J is split at `sbpR` into J0 J1.
        Then we have

            I0 I1 S J0 J1

        and the tandem duplication produces

            I0 I1 S J0 I1 S J0 J1
        """
        self.afterL = copy.deepcopy(self.beforeL)

        lenI1 = self.beforeL.sc2 - self.sbpL #number of (intergene) muclotides
        lenJ0 = self.sbpR - self.beforeR.sc1 #number of (intergene) muclotides
        lenJ1 = self.beforeR.sc2 - self.sbpR #number of (intergene) muclotides
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
    
        position = self.beforeL.position + self.numintergenes + 1
        self.afterC = Interval(tcenterstart, tcenterend,
                               scenterstart, scenterend, position, 'I',
                               tcenterbreak, scenterbreak)
        position += self.numintergenes + 1
        self.afterR = Interval(trightstart, trightend,
                               srightstart, srightend, position, 'I')

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
            return self.assertS(self.beforeL.sc1 + (sc - self.afterL.sc1))

        elif self.afterC.inSpecific(sc):
            if sc <= self.afterC.s_bp:                  # J0
                return self.assertS(self.beforeR.sc1 + (sc - self.afterC.sc1))
            else:                                       # I1
                return self.assertS(self.beforeL.sc2 - (self.afterC.sc2 - sc))

        elif self.afterR.inSpecific(sc):                # J0 J1
            return self.assertS(self.beforeR.sc1 + (sc - self.afterR.sc1))

        elif sc > self.afterC.sc2:                      # to the right
            offset = sc - self.afterC.s_bp
            if self.wraps() and offset > self.swraplen - self.sbpL:
                return self.assertS(sc - (self.lenSs - 1) - self.afterC.specificLen())
            else:
                return self.assertS(self.sbpL + offset)
        else:                                           # to the left
            return self.assertS(sc)

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
            return self.assertT(self.beforeL.tc1 + (tc - self.afterL.tc1))

        elif self.afterC.inTotal(tc):
            if tc <= self.afterC.t_bp:          # J0
                return self.assertT(self.beforeR.tc1 + (tc - self.afterC.tc1))
            else:                               # I1
                return self.assertT(self.beforeL.tc2 - (self.afterC.tc2 - tc))

        elif self.afterR.inTotal(tc):           # J0 J1
            return self.assertT(self.beforeR.tc1 + (tc - self.afterR.tc1))

        elif tc > self.afterC.tc2:              # to the right
            offset = tc - self.afterC.t_bp
            if self.wraps() and offset > self.twraplen - self.tbpL:
                return self.assertT(tc - (self.lenSt - 1) - self.afterC.specificLen())
            else:
                return self.assertT(self.tbpL + offset)
        else:                                   # to the left
            return self.assertT(tc)

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
               f'{self.swraplen} {self.lineage} {self.time}'
