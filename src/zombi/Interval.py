from typing import Tuple

from . import T_PAIR

class Interval:
    """
    The breakpoint interval of a gene or an intergene.
    A Gene or Intergene interval holding both total and specific coordinates, 
    along with possibly a breakpoint coordinate contained inside it. Total
    coordinates are with respect to all nucleotides (Gene or Intergene), whereas
    specific coordinates consider only the nucleotides from `itype`, ignoring
    the others.

    Attributes
    ----------
    tc1: int
        left total breakpoint coordinate of interval
    tc2: int
        right total breakpoint coordinate of interval, inclusive
    sc1: int
        left intergene specific breakpoint coordinate of interval
    sc2: int
        right intergene specific breakpoint coordinate of interval, inclusive
    position: int
        the position of the Gene or Intergene in its list
    itype: str
        one of 'G' or 'I' for Gene or Intergene
    t_bp: int
        the total breakpoint coordinate
    s_bp: int
        the specific breakpoint coordinate
    """
    def __init__(self, tc1: int, tc2: int, sc1: int, sc2: int, index: int,
                 itype: str, total=-1, specific=-1):
        """
        Create a new Interval that may or may not contain breakpoint
        coordinates.

        Parameters
        ----------
        tc1 : int
            first total coordinate
        tc2 : int
            second total coordinate, inclusive
        spc1 : int
            first specific coordinate
        spc2 : int
            second specific coordinate, inclusive
        index : int
            the position of the Interval in its respective list
        itype : str
            one of 'I' or 'G'
        total : int, optional
            total coordinate in the interval corresponding to a breakpoint, by
            default -1
        specific : int, optional
            specific coordinate in the interval corresponding to `total`, by
            default -1
        """
        self.tc1: int = tc1
        self.tc2: int = tc2
        self.sc1: int = sc1
        self.sc2: int = sc2
        self.position: int = index
        self.itype: str = itype
        self.t_bp: int = total
        self.s_bp: int = specific

            #Sanity checks:
        assert tc2 - tc1 == sc2 - sc1, f'{tc2-tc1} != {sc2-sc1}'
        if total >= 0:
            assert tc1 <= total <= tc2, 'total breakpoint outside of range'
        if specific >= 0:
            assert sc1 <= specific <= sc2, 'specific breakpoint outside of range'
            if total >= 0:
                assert tc2 - total == sc2 - specific, 'breakpoint mismatch'

    def asTuple(self) -> Tuple[int, int, int, int, int, str]:
        """
        Convert this Location into a tuple.

        Returns
        -------
        Tuple[int, int, int, int, int, str]
            (tc1, tc2, sc1, sc2, position, itype)
        """
        return self.tc1, self.tc2, self.sc1, self.sc2, self.position, self.itype

    def totalPair(self) -> T_PAIR:
        return self.tc1, self.tc2

    def specificPair(self) -> T_PAIR:
        return self.sc1, self.sc2

    def specificToTotal(self, sc: int) -> int:
        return self.tc1 + (sc - self.sc1)

    def totalToSpecific(self, tc: int) -> int:
        return self.sc1 + (tc - self.tc1)

    def inTotal(self, t_coord:int) -> bool:
        """
        Is the given total coordinate inside this Interval?
        """
        return self.tc1 <= t_coord <= self.tc2

    def inSpecific(self, s_coord:int) -> bool:
        """
        Is the given specific coordinate inside this Interval?
        """
        return self.sc1 <= s_coord <= self.sc2

    def specificLen(self) -> int:
        """
        Return the number of specific breakpoint coordinates in this interval.
        """
        return self.sc2 - self.sc1 + 1

    def totalLen(self) -> int:
        """
        Return the number of total breakpoint coordinates in this interval.
        """
        return self.tc2 - self.tc1 + 1

    def isIntergenic(self) -> bool:
        return self.itype == 'I'

    def splitSpecific(self) -> Tuple[T_PAIR, T_PAIR]:
        """
        Get the halves of the specific interval split by `s_bp`.
        """
        return (self.sc1, self.s_bp), (self.s_bp, self.sc2)

    def splitTotal(self) -> Tuple[T_PAIR, T_PAIR]:
        """
        Get the halves of the specific interval split by `t_bp`.
        """
        return (self.tc1, self.t_bp), (self.t_bp, self.tc2)

    def __len__(self):
        return self.sc2 - self.sc1

    def __str__(self):
        return f'{self.tc1} {self.tc2} {self.sc1} {self.sc2} {self.position} {self.itype}'

    def __repr__(self):
        return f'{self.tc1} {self.tc2} {self.sc1} {self.sc2} {self.s_bp}'
