from typing import Tuple

from . import T_PAIR

class Interval:
    """
    A Gene or Intergene interval holding both total and specific coordinates, 
    along with possibly a breakpoint coordinate contained inside it. Total
    coordinates are with repect to all nucleotides (Gene or Intergene), whereas
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
    t_breakpoint: int
        the total breakpoint coordinate
    s_breakpoint: int
        the specific breakpoint coordinate
    """
    def __init__(self, tc1: int, tc2: int, sc1: int, sc2: int, index: int,
                 itype: str, total=0, specific=0):
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
            default 0
        specific : int, optional
            specific coordinate in the interval corresponding to `total`, by default 0
        """
        self.tc1 = tc1
        self.tc2 = tc2
        self.sc1 = sc1
        self.sc2 = sc2
        self.position = index
        self.itype = itype
        self.t_breakpoint = total
        self.s_breakpoint = specific

            #Sanity checks:
        assert tc2 - tc1 == sc2 - sc1, print(f'{tc2-tc1} != {sc2-sc1}')
        if total:
            assert tc1 <= total <= tc2
        if specific:
            assert sc1 <= specific <= sc2
            if total:
                assert tc2 - total == sc2 - specific

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

    def containsTotal(self, t_coord:int) -> bool:
        return self.tc1 <= t_coord <= self.tc2

    def containsSpecific(self, s_coord:int) -> bool:
        return self.sc1 <= s_coord <= self.sc2

    def isIntergenic(self) -> bool:
        return self.itype == 'I'

    def splitSpecific(self) -> Tuple[T_PAIR, T_PAIR]:
        """
        Get the halves of the specific interval split by `s_breakpoint`.
        """
        return (self.sc1, self.s_breakpoint), (self.s_breakpoint, self.sc2)

    def splitTotal(self) -> Tuple[T_PAIR, T_PAIR]:
        """
        Get the halves of the specific interval split by `t_breakpoint`.
        """
        return (self.tc1, self.t_breakpoint), (self.t_breakpoint, self.tc2)

    def __len__(self):
        return self.sc2 - self.sc1

    def __str__(self):
        return f'{self.tc1} {self.tc2} {self.sc1} {self.sc2} {self.position} {self.itype}'

    def __repr__(self):
        return str(self)