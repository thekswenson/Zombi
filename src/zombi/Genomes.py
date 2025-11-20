"""
Classes to represent the associated parts of a genome.  Genomes have
Chromosomes, which contain Genes and Intergenes.  Intergenes are composed of
Divisions, which are grouped into DivisionFamilies.  Genes are grouped into
GeneFamilies.

Each GeneFamily and DivisionFamily keeps track of its history of events on its
own "gene tree" (these trees are implicitly held in events lists, rather than
explicitly being represented).

A GeneFamily or DivisionFamily has nodes names like "n16;2" in the files (or
"n16_2" in the code), where "n16" is the name of the species tree node where the
gene/division resides, and "2" is a unique identifier for the node in that
gene/division tree (i.e. there exactly one occurrence of "2" in the family).
"""
import abc
import itertools
import sys
import ete3
import networkx as nx
import numpy as np

from functools import reduce
from typing import Iterable, Union, cast
from enum import Enum, auto


from . import AuxiliarFunctions as af
from . import ReconciledTree as RT
from . import T_PAIR
from .Interval import ITYPE, Interval
from .Events import BACKWARDS, FORWARDS, GenomeCoordEvent
from .Events import LOSS, ORIG, POS, TDUP, DUP, FER, INV
from .Random import G_RNG, G_NPRNG


# Directions:
class T_DIR(Enum):
    """ The direction to extend from an index. """
    RIGHT = auto()
    LEFT = auto()

    
#T_DIR = bool #: the direction
#
#RIGHT = False
#LEFT = True



class GeneFamily:
    """
    Represent a gene family which knows its history of events on a gene tree.

    Attributes
    ----------
    events: list[tuple[str, str, str]]
        list of events (time, type, location) where time is a float , type is
        the single capital letter representing the event type (e.g. 'O', 'I',
        etc.), and location is a ';' delimited list of node names with gene ids
        (e.g. n16;2;n8;35).  See the Zombi wiki for details on the location.
    """
    def __init__(self, identifier, time):

        self.identifier = identifier    #unique integer
        self.origin = time

        self.gff_id = ''                #unique ID from the gff file

        self.genes: list[Gene] = []
        self.events: list[tuple[float, str, str]] = []
        self.event_counter = 0  # Each time that the family is modified in any form, we have to update the event counter
        self.gene_ids_counter = 0

        self.length = 0

        self.rates = dict()           # Only in Gm mode
        self.initial_orientation = "" # The initial orientation of the family


    def append_event(self, time: float, event: str, genes: str):

        self.events.append((time, event, genes))


    def generate_tree(self) -> tuple[str, str, str]:

        def find_descendant(surviving_nodes, node):

            found = 0
            mynode = surviving_nodes[node]["descendant"]

            while found == 0:

                if surviving_nodes[mynode]["state"] == 1:
                    found = 1
                else:
                    mynode = surviving_nodes[mynode]["descendant"]

            return mynode

        # Eric's algorithm

        # First we will iterate the events from the end

        events = self.events

        surviving_nodes = dict()
        times = dict()

        family_size = 0

        for current_time, event, nodes in reversed(events):

            if event == "F":

                nodename = nodes.replace(";","_")
                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 1, "descendant": "None"}

                family_size += 1

            elif event == "E" or event == LOSS:

                nodename = nodes.replace(";", "_")

                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 0, "descendant": "None"}

            elif event == "S" or event == TDUP or event == DUP or event == FER:

                p, g0, c1, g1, c2, g2 = nodes.split(";")

                pnodename = p + "_" + g0
                c1nodename = c1 + "_" + g1
                c2nodename = c2 + "_" + g2

                times[pnodename] = float(current_time)

                if surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == 1:

                    surviving_nodes[pnodename] = {"state": 1, "descendant": c1nodename + ";" + c2nodename}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == 0:

                    surviving_nodes[pnodename] = {"state": 0, "descendant": "None"}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode1 = find_descendant(surviving_nodes, c1nodename)
                    mynode2 = find_descendant(surviving_nodes, c2nodename)

                    surviving_nodes[pnodename] = {"state": 1, "descendant": mynode1 + ";" + mynode2}

                elif surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == 0:

                    surviving_nodes[pnodename] = {"state": -1, "descendant": c1nodename}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == 1:

                    surviving_nodes[pnodename] = {"state": -1, "descendant": c2nodename}

                elif surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2nodename)
                    surviving_nodes[pnodename] = {"state": 1, "descendant": c1nodename + ";" + mynode}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == 1:

                    mynode = find_descendant(surviving_nodes, c1nodename)
                    surviving_nodes[pnodename] = {"state": 1, "descendant": mynode + ";" + c2nodename}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == 0:

                    mynode = find_descendant(surviving_nodes, c1nodename)
                    surviving_nodes[pnodename] = {"state": -1, "descendant": mynode}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2nodename)
                    surviving_nodes[pnodename] = {"state": -1, "descendant": mynode}
                

        extanttree = RT.ReconciledTree()
        completetree = RT.ReconciledTree()

        eroot = extanttree.get_tree_root()
        eroot.name = ""

        wquick_nodes = dict()
        equick_nodes = dict()

        for values in events:

            current_time, event, nodes = values

            if event == ORIG:

                wroot = completetree.get_tree_root()
                wroot.name = nodes + "_1"
                wquick_nodes[wroot.name] = wroot

            if event == LOSS or event == "E":

                p, g0 = nodes.split(";")
                pnodename = p + "_" + g0
                mynode = wquick_nodes[pnodename]
                e = RT.RecEvent(LOSS, p, int(float(current_time)))
                mynode.addEvent(e, append=True)

            if event == "F":

                p, g0 = nodes.split(";")
                pnodename = p + "_" + g0
                mynode = wquick_nodes[pnodename]
                e = RT.RecEvent(POS, p, int(float(current_time)))
                mynode.addEvent(e, append=True)

            if event == "S" or event == TDUP or event == DUP or event == FER:

                p, g0, c1, g1, c2, g2 = nodes.split(";")
                pnodename = p + "_" + g0
                c1nodename = c1 + "_" + g1
                c2nodename = c2 + "_" + g2

                mynode = wquick_nodes[pnodename]
                myc1 = mynode.add_child()
                myc2 = mynode.add_child()
                myc1.name = c1nodename
                myc2.name = c2nodename
                myc1.dist = times[c1nodename] - times[pnodename]
                myc2.dist = times[c2nodename] - times[pnodename]

                wquick_nodes[c1nodename] = myc1
                wquick_nodes[c2nodename] = myc2

                state = surviving_nodes[pnodename]["state"]

                ### Now we add the reconciled events

                e = RT.RecEvent(event, p, int(float(current_time)))
                mynode.addEvent(e, append=True)

                if state == 1:  # Now the extant tree

                    c1name, c2name = surviving_nodes[pnodename]["descendant"].split(";")

                    if eroot.name == "":
                        eroot.name = pnodename
                        equick_nodes[pnodename] = eroot

                    mynode = equick_nodes[pnodename]

                    myc1 = mynode.add_child()
                    myc2 = mynode.add_child()

                    myc1.name = c1name
                    myc2.name = c2name

                    myc1.dist = times[c1name] - times[pnodename]
                    myc2.dist = times[c2name] - times[pnodename]

                    equick_nodes[c1name] = myc1
                    equick_nodes[c2name] = myc2

        if family_size == 0:
            extanttreestr = ";"

        elif family_size == 1:
            extanttreestr = [k for k, v in surviving_nodes.items() if v["state"] == 1 and v["descendant"] == "None"][
                             0] + ";"
        else:
            extanttreestr = extanttree.write(format=1, format_root_node=True)


        rec = completetree.getTreeRecPhyloXML()

        if len(completetree) == 0:
            completetreestr = ";"
        elif len(completetree) == 1:
            completetreestr = completetree.get_leaves()[0].name + ";"
        else:
            completetreestr = completetree.write(format=1, format_root_node=True)

        assert isinstance(completetreestr, str)
        assert isinstance(extanttreestr, str)
        return completetreestr, extanttreestr, rec


    def generate_oldtree(self):

        tree = ete3.Tree()

        current_time, event, nodes = self.events[0]

        sp = tree.get_tree_root()
        sp.name = nodes + "_1"
        sp.add_feature("is_active", True)

        elapsed_time = float(current_time)

        for current_time, event, nodes in self.events[1:]:

            elapsed_time = float(current_time) - elapsed_time
            active_nodes = [x for x in tree.get_leaves() if x.is_active == True]    #type: ignore
            for node in active_nodes:
                node.dist += elapsed_time
            elapsed_time = float(current_time)

            if event == "S":

                sp, gp, c1, g1, c2, g2 = nodes.split(";")

                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False        #type: ignore

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == "E":
                sp, gp = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False        #type: ignore

            elif event == LOSS:
                sp, gp = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False        #type: ignore

            elif event == TDUP or event == DUP:

                sp, gp, c1, g1, c2, g2 = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname

                mynode.is_active = False        #type: ignore

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == FER:
                sp, gp, c1, g1, c2, g2 = nodes.split(";")

                myname = sp + "_" + gp

                mynode = tree & myname
                mynode.is_active = False        #type: ignore

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == "F":
                break

        complete_tree = tree.write(format=1, format_root_node=True)
        active_nodes = [x for x in tree.get_leaves() if x.is_active == True]    #type: ignore

        if len(active_nodes) < 3:
            pruned_tree = None

        else:
            tree.prune(active_nodes, preserve_branch_length=True)
            pruned_tree = tree.write(format=1, format_root_node=True)

        return complete_tree, pruned_tree

    def obtain_new_gene_id(self):
        self.gene_ids_counter += 1
        return self.gene_ids_counter

    def __str__(self):

        return "GeneFamily_" + str(self.identifier) + ";" + ";".join([str(x) for x in self.genes])

        #return ";".join(map(str, self.genes))

    def __len__(self):

        return len([x for x in self.genes if x.active == True])

    def __iter__(self):
        for gene in self.genes:
            yield gene




class Gene:
    """
    Attributes
    ----------
    family: str
        the name of the gene family
    gene_id: int
        the unique identifier of the gene within the family
    length: int
        the length of the gene
    start: int
        the starting index of the gene (in gene order), python indexed
    end: int
        the ending index of the gene, python indexed (non-inclusive)
    orientation: str
        the orientation of the gene ("+" or "-")
    total_flanking: T_PAIR
        these are the coordinates of the breakpoints at either end of the gene.
        (i, j) where i is the total number of breakpoints before this gene
        (assuming no intergene before the first gene), and j is i plus the
        length of this gene
    specific_flanking: T_PAIR
        these are the coordinates of the breakpoints at either end of the gene,
        while only counting breakpoints that touch a gene.
        (i, j) where i is the number of gene breakpoints (not intergene) before
        this gene excluding all the breakpoints from the intergenes, and j is i
        plus the length of this gene
    """

    def __init__(self):

        self.active = True    # Inactive when no longer on the current branch,
                              # or replaced (warning even events like transfer,
                              # duplication, speciation will deactivate a gene)
        self.orientation = ""
        self.family = ""      # FIX this variable should have the same name that the division one
        self.gene_id = -1     # FIX this variable should have the same name that the division one
        self.sequence = ""
        self.species = ""     #: Also known as "lineage" in some places
        self.importance = 0
        self.length = 0
        self.start: int              #: pythonic (inclusive start, 0 indexed)
        self.end: int                #: pythonic (non-inclusive end)
        self.total_flanking: T_PAIR         #: not pythonic (both inclusive)
        self.specific_flanking: T_PAIR      #: not pythonic (both inclusive)
        self.ptype = "Gene"                 # For debugging purposes

    def determine_orientation(self):

        if G_NPRNG().binomial(1,0.5):
            self.orientation = "+"

        else:
            self.orientation = "-"

    def change_orientation(self):

        if self.orientation == "+":
            self.orientation = "-"
        elif self.orientation == "-":
            self.orientation = "+"
        else:
            raise(Exception('unexpected orientation "{self.orientation}"'))

    def __str__(self):

        myname = "_".join(map(str, (self.family, self.length)))
        #myname = "_".join(map(str, (self.gene_family, self.orientation)))
        #myname = str(self.gene_family) + "_" + str(self.gene_id)
        #myname = "_".join(map(str, (self.gene_family, self.length)))
        return myname

    def __repr__(self):
        return f'{self.family}_{self.gene_id}'

    def __len__(self):
        return self.length


class Division:
    """
    A class that represents all the subdivisions of an intergene
    
    Attributes
    ----------
    events: list
        a list with all the events affecting this division
    """
    def __init__(self, identity: int, division_family: int, specific_flanking: T_PAIR|None=None):
        
        self.orientation = "+"
        self.identity = identity # Specific identifier of this division
        self.family = division_family # The name of the division family
        self.specific_flanking: T_PAIR
        self.length = 0
        if specific_flanking:
            self.specific_flanking = specific_flanking   #: not pythonic (both inclusive)
            self.set_length()
        self.ptype = "Divi" # Piece type, for debugging purposes
        self.species = ""
        self.total_flanking: T_PAIR

        self.initial_sequence = "" # The initial sequence if the gene comes from a pseudogenization

    def change_orientation(self):

        if self.orientation == "+":
            self.orientation = "-"
        elif self.orientation == "-":
            self.orientation = "+"

    def set_length(self):
        self.length = int(abs(self.specific_flanking[1] - self.specific_flanking[0]))
        
    def __str__(self):
        return  str(self.family) + "_" + str(self.length) + "_" + str(self.specific_flanking)

    def __len__(self) -> int:
        return self.length


class Intergene:
    """
    Attributes
    ----------
    length: int
        the length of the intergene (in nucleotides, not breakpoints)
    total_flanking: T_PAIR
        these are the coordinates of the breakpoints at either end of the
        intergene.  (i, j) where i is the total number of breakpoints before
        this intergene, and j is i plus the length of this intergene
    specific_flanking: T_PAIR
        these are the coordinates of the breakpoints at either end of the
        intergene, while only counting breakpoints that touch an intergene.
        (i, j) where i is the number of intergenic breakpoints before this
        intergene, and j is i plus the length of this intergene
    """

    def __init__(self, length = 0):
        assert length >= 0, f"Intergene length {length} must be non-negative!"
        if length:
            self.length = length
        else:
            self.length = 0

        self.total_flanking: T_PAIR             #: not pythonic (both inclusive)
        self.specific_flanking: T_PAIR          #: not pythonic (both inclusive)
        self.id = 0                             # Only for debugging purposes
        self.divisions: list[Division] = list() # List containing the divisions

    @property
    def sc1(self):
        return self.specific_flanking[0]

    @property
    def sc2(self):
        return self.specific_flanking[1]

    @property
    def tc1(self):
        return self.total_flanking[0]

    @property
    def tc2(self):
        return self.total_flanking[1]

    def inTotal(self, tc:int) -> bool:
        return self.total_flanking[0] <= tc <= self.total_flanking[1]

    def inSpecific(self, sc:int) -> bool:
        return self.specific_flanking[0] <= sc <= self.specific_flanking[1]

    def create_division(self, identity: int, division_family, specific_flanking):
        """
        """
        division = Division(identity, division_family, specific_flanking)
        self.divisions.append(division)
        return division
        
    def __str__(self):

        #return "(" + str(self.length) + ")"
        #return "I_" + str(self.length)
        #return "I_" + str(self.id) + "_" + str(self.length)
        #if not self.specific_flanking:
        #    raise(Exception("Specific flanking indices not innitilized yet!"))
        #return "I_" + str(self.specific_flanking) 
        if not self.total_flanking:
            raise(Exception("Total flanking indices not innitilized yet!"))
        return "I_" + str(self.total_flanking) 
    
    def __iter__(self):

        for x in self.divisions:
            yield x
   
    def __len__(self):
       return self.length


class DivisionFamily:
    """
    A division family is akin to a GeneFamily. It stores all the events in an
    ensemble of homologous divisions.

    Attributes
    ----------
    events: list[tuple[str, str, str]]
        list of events (time, type, location) where time is a float , type is
        the single capital letter representing the event type (e.g. 'O', 'E',
        'L', 'F', 'L', 'D', 'U', 'T', 'S'), and location is a ';' delimited list
        of node names with division ids (e.g. n16;2;n8;35).  See the Zombi wiki
        for details on the location.
    """
    def __init__(self, id: int, initial_specific_flanking: T_PAIR):

        self.id = id
        self.events: list[tuple[float, str, str]] = list()
        self.gene_ids_counter = 1
        self.initial_flanking  = initial_specific_flanking
        self.initial_orientation = ""

    def append_event(self, time: float, event: str, divisions: str):
        self.events.append((time, event, divisions))
    
    def obtain_new_identifier(self) -> int:
        self.gene_ids_counter += 1
        return self.gene_ids_counter


    def generate_tree(self) -> tuple[str, str, str]:

        def find_descendant(surviving_nodes, node):

            found = False
            mynode = surviving_nodes[node]["descendant"]

            while not found:

                if surviving_nodes[mynode]["state"] == 1:
                    found = True
                else:
                    mynode = surviving_nodes[mynode]["descendant"]

            return mynode

        # Eric's algorithm

        # First we will iterate the events from the end

        events = self.events

        surviving_nodes = dict()
        times = dict()

        family_size = 0

        for current_time, event, nodes in reversed(events):

            if event == "F":

                nodename = nodes.replace(";","_")
                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 1, "descendant": "None"}

                family_size += 1

            elif event == "E" or event == LOSS:

                nodename = nodes.replace(";", "_")

                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 0, "descendant": "None"}

            elif event == "S" or event == TDUP or event == DUP or event == FER:

                p, g0, c1, g1, c2, g2 = nodes.split(";")

                pnodename = p + "_" + g0
                c1nodename = c1 + "_" + g1
                c2nodename = c2 + "_" + g2

                times[pnodename] = float(current_time)

                if surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == 1:

                    surviving_nodes[pnodename] = {"state": 1, "descendant": c1nodename + ";" + c2nodename}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == 0:

                    surviving_nodes[pnodename] = {"state": 0, "descendant": "None"}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode1 = find_descendant(surviving_nodes, c1nodename)
                    mynode2 = find_descendant(surviving_nodes, c2nodename)

                    surviving_nodes[pnodename] = {"state": 1, "descendant": mynode1 + ";" + mynode2}

                elif surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == 0:

                    surviving_nodes[pnodename] = {"state": -1, "descendant": c1nodename}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == 1:

                    surviving_nodes[pnodename] = {"state": -1, "descendant": c2nodename}

                elif surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2nodename)
                    surviving_nodes[pnodename] = {"state": 1, "descendant": c1nodename + ";" + mynode}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == 1:

                    mynode = find_descendant(surviving_nodes, c1nodename)
                    surviving_nodes[pnodename] = {"state": 1, "descendant": mynode + ";" + c2nodename}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == 0:

                    mynode = find_descendant(surviving_nodes, c1nodename)
                    surviving_nodes[pnodename] = {"state": -1, "descendant": mynode}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2nodename)
                    surviving_nodes[pnodename] = {"state": -1, "descendant": mynode}

            elif event == ORIG:
                pass

            else:
                raise(Exception(f"Unknown event type {event}!"))

        extanttree = RT.ReconciledTree()
        completetree = RT.ReconciledTree()

        eroot = extanttree.get_tree_root()
        eroot.name = ""

        wquick_nodes = dict()
        equick_nodes = dict()

        for values in events:

            current_time, event, nodes = values

            if event == ORIG:

                wroot = completetree.get_tree_root()
                wroot.name = nodes + "_1"
                wquick_nodes[wroot.name] = wroot

            if event == LOSS or event == "E":

                p, g0 = nodes.split(";")
                pnodename = p + "_" + g0
                mynode = wquick_nodes[pnodename]
                e = RT.RecEvent(LOSS, p, int(float(current_time)))
                mynode.addEvent(e, append=True)

            if event == "F":

                p, g0 = nodes.split(";")
                pnodename = p + "_" + g0
                mynode = wquick_nodes[pnodename]
                e = RT.RecEvent(POS, p, int(float(current_time)))
                mynode.addEvent(e, append=True)

            if event == "S" or event == TDUP or event == DUP or event == FER:

                p, g0, c1, g1, c2, g2 = nodes.split(";")
                pnodename = p + "_" + g0
                c1nodename = c1 + "_" + g1
                c2nodename = c2 + "_" + g2

                mynode = wquick_nodes[pnodename]
                myc1 = mynode.add_child()
                myc2 = mynode.add_child()
                myc1.name = c1nodename
                myc2.name = c2nodename
                myc1.dist = times[c1nodename] - times[pnodename]
                myc2.dist = times[c2nodename] - times[pnodename]

                wquick_nodes[c1nodename] = myc1
                wquick_nodes[c2nodename] = myc2

                state = surviving_nodes[pnodename]["state"]

                ### Now we add the reconciled events

                e = RT.RecEvent(event, p, int(float(current_time)))
                mynode.addEvent(e, append=True)

                if state == 1:  # Now the extant tree

                    c1name, c2name = surviving_nodes[pnodename]["descendant"].split(";")

                    if eroot.name == "":
                        eroot.name = pnodename
                        equick_nodes[pnodename] = eroot

                    mynode = equick_nodes[pnodename]

                    myc1 = mynode.add_child()
                    myc2 = mynode.add_child()

                    myc1.name = c1name
                    myc2.name = c2name

                    myc1.dist = times[c1name] - times[pnodename]
                    myc2.dist = times[c2name] - times[pnodename]

                    equick_nodes[c1name] = myc1
                    equick_nodes[c2name] = myc2

        if family_size == 0:
            extanttreestr = ";"

        elif family_size == 1:
            extanttreestr = [k for k, v in surviving_nodes.items()
                             if v["state"] == 1 and v["descendant"] == "None"][0] + ";"
        else:
            extanttreestr = extanttree.write(format=1, format_root_node=True)

        rec = completetree.getTreeRecPhyloXML()

        if len(completetree) == 0:
            completetreestr = ";"
        elif len(completetree) == 1:
            completetreestr = completetree.get_leaves()[0].name + ";"
        else:
            completetreestr = completetree.write(format=1, format_root_node=True)

        assert isinstance(completetreestr, str)
        assert isinstance(extanttreestr, str)
        return completetreestr, extanttreestr, rec
        

    def __repr__(self):
        return str(self.initial_flanking)

    def __len__(self):
        return int(self.initial_flanking[1] -  self.initial_flanking[0])

    
  

class Chromosome(abc.ABC):
    """
    A chromosome that knows its genes and intergenes, as well as its
    `map_of_locations`, which is the representation of the chromosome as a
    list of Intervals representing the (interleaved) Genes and Intergenes.

    Attributes
    ----------
    map_of_locations: List[Location]
        list of Locations, which represent gene or intergene regions
    _num_nucleotides: int
        the length of the chromosome (in nucleotides)
    genes: list[Gene]
        ordered list of genes in this chromosome
    intergenes: list[Intergene]
        ordered list of intergenes in this chromosome
    shape: str
        one of "L" or "C" for linear or circular
    event_history: GenomeCoordEvent
        list of genome events (e.g. INV, TDUP, etc.) that have modified this
        chromosome, with the genome coordinates of where they occurred
    geneorder_history: GenomeEvent
        list of genome events (e.g. INV, TDUP, etc.) that have modified this
        chromosome, with the breakpoint (gene-order) positions of where they
        occurred
    pieces: list[Gene|Division]
        In the F mode, an ordered list of genes and divisions that is used to
        reconstruct the rearrangement sequence once all divisions are known.
    """
    def __init__(self, name: str, num_nucleotides=0):

        self.name = name
        self.has_intergenes = False
        self.intergenes: list[Intergene] = list()
        self.genes: list[Gene] = list()
        self.shape = ""
        self._num_nucleotides = num_nucleotides  # length in nucleotides

        self.map_of_locations: list[Interval] = []

        self.total_rates = 0

        self.event_history: list[GenomeCoordEvent] = []
        #self.geneorder_history: list[GeneOrderEvent] = []

        self.pieces: list[Gene|Division] = []  # In the F mode, keeps a list of genes and divisions


    def obtain_total_itergenic_length(self):

        total_length = 0
        for intergene in self.intergenes:
            total_length += intergene.length
        return total_length


    def select_random_position(self, blacklist: list[int]=[]) -> int:
        """
        Get a random gene position.

        Parameters
        ----------
        blacklist: list[int]
            A list of positions that are not allowed to be selected.
        """
        if blacklist:
            positions = [i for i in range(len(self.genes)) if i not in blacklist]
            if not positions:
                raise ValueError("All positions blacklisted!")
            return G_NPRNG().choice(sorted(positions))

        return int(G_NPRNG().integers(len(self.genes)))


    def fill_pieces(self):
        """
        Once we know the divisions and genes, we can use this function
        to fill them into the list of pieces in the right order,
        and obtain the total flankings of all pieces
        """
        for gene, intergene in zip(self.genes, self.intergenes):
            self.pieces.append(gene)
            right_bp = gene.total_flanking[1] 
            for division in intergene:
                division.total_flanking = (right_bp, right_bp + len(division))
                right_bp =  right_bp + len(division)
                self.pieces.append(division) 

    def get_index_gene(self, search_gene):
        """
        Obtain the index of the gene in the pieces
        The index refers only at the genes, so
        intergenic divisions are ignored
        """
        genes = [piece for piece in self.pieces if isinstance(piece, Gene)]
        return (genes.index(search_gene))

    def print_pieces(self):       
        """
        For debugging purposes. Print all pieces in the genome
        """
        for piece in self.pieces:
            if isinstance(piece, Gene):
                print(piece.ptype, piece.total_flanking, piece.length, piece.orientation, piece.family)
            else:
                assert isinstance(piece, Division)
                print(piece.ptype, piece.total_flanking, piece.length, piece.orientation, piece.family, piece.specific_flanking, piece.initial_sequence) 


    def update_coordinates(self):

        """
        Update the coordinates of the genes and the divisions
        """
        right_bp = len(self.pieces[0])
        self.pieces[0].total_flanking = (0, right_bp)
        for piece in self.pieces[1:]:
            piece.total_flanking = (right_bp, right_bp + len(piece))
            right_bp =  right_bp + len(piece)

        self.update_specific_coordinates()


    def update_specific_coordinates(self):
        """
        Get the specific coordinates of the intergenes
        """
        first_intergene = True
        adjacent_gene = False
        adjacent_division = False
        coordinate_adjustment = 0 

        self.specific2total = dict()

        right_bp = 0    #for pylance
        for piece in self.pieces:
    
            if isinstance(piece, Gene):
                adjacent_division = False
                if adjacent_gene == True:
                    coordinate_adjustment += 1
                adjacent_gene = True
                continue

            elif isinstance(piece, Division):
                adjacent_gene = False

                if first_intergene == True:
                    
                    adjacent_division = True
                    first_intergene = False
                    piece.specific_flanking = (coordinate_adjustment, coordinate_adjustment + len(piece))
                    right_bp = coordinate_adjustment + len(piece) 
                    coordinate_adjustment = 0
                    
                else:

                    if adjacent_division == False:
                        right_bp += 1
                    adjacent_division = True
                    piece.specific_flanking = (right_bp +  coordinate_adjustment, right_bp + coordinate_adjustment + len(piece))
                    right_bp = right_bp + coordinate_adjustment + len(piece)
                    coordinate_adjustment = 0

                assert isinstance(piece.total_flanking, tuple), "Total flanking not set yet!"
                self.specific2total[piece.specific_flanking[0]] = piece.total_flanking[0]
                self.specific2total[piece.specific_flanking[1]] = piece.total_flanking[1]
        
        

    def update_flankings(self):
        """
        Set the "flanking" breakpoint intervals for each of the genes and
        intergenes based on their lengths.
        """

        if self.has_intergenes:

            self.genes[0].total_flanking = (0, self.genes[0].length)
            self.genes[0].specific_flanking = (0, self.genes[0].length)
            self.intergenes[0].total_flanking = (self.genes[0].total_flanking[1],
                                                 self.genes[0].total_flanking[1] + self.intergenes[0].length)
            self.intergenes[0].specific_flanking = (0, self.intergenes[0].length)

            for i in range(len(self.genes)):

                if i == 0:
                    continue

                lb = self.intergenes[i-1].total_flanking[1]
                ub = lb + self.genes[i].length

                lbg = self.genes[i-1].specific_flanking[1] + 1
                ubg = lbg + self.genes[i].length

                self.genes[i].total_flanking = (lb, ub)
                self.genes[i].specific_flanking = (lbg, ubg)

                lbi = self.intergenes[i-1].specific_flanking[1] + 1
                ubi = lbi + self.intergenes[i].length

                self.intergenes[i].total_flanking = (ub, ub + self.intergenes[i].length)
                self.intergenes[i].specific_flanking = (lbi, ubi)

            #self.intergenes[i].total_flanking = (ub, 0)
            

    def update_locations(self):
        """
        Setup the `map_of_locations` list which will contain all of the genes
        and intergenes interleaved, in the order in which they appear in the
        genome. In the process, set the "flanking" endpoints for each of the
        genes and intergenes (based on their length).
        """
        self.update_flankings()
        self.map_of_locations = list()

        for i in range(len(self.genes)):

            tc1 = self.genes[i].total_flanking[0]
            tc2 = self.genes[i].total_flanking[1]
            sc1 = self.genes[i].specific_flanking[0]
            sc2 = self.genes[i].specific_flanking[1]
            self.map_of_locations.append(Interval(tc1, tc2, sc1, sc2, i,
                                                  ITYPE.GENE, gene=self.genes[i]))

            tc1 = self.intergenes[i].total_flanking[0]
            tc2 = self.intergenes[i].total_flanking[1]
            sc1 = self.intergenes[i].specific_flanking[0]
            sc2 = self.intergenes[i].specific_flanking[1]
            self.map_of_locations.append(Interval(tc1, tc2, sc1, sc2, i,
                                                  ITYPE.INTERGENE))
    
    def update_flankings_divisions(self):

        """
        After an event, updates the coordinates of the divisions. 
        """

        for intergene in self.iter_intergenes():

            if len(intergene.divisions) == 0:
                continue

            lb, rb = intergene.specific_flanking # the specific flankings of the intergene
            first_division = intergene.divisions[0]
            
            first_division.specific_flanking = (lb, lb + len(first_division))
            
            for i, division in enumerate(intergene.divisions):
                if i == 0:
                    continue
                previous_division = intergene.divisions[i-1]
                division.specific_flanking = (previous_division.specific_flanking[1], previous_division.specific_flanking[1] + len(division))


    def select_random_coordinate_in_intergenic_regions(self, exclude: list[int]=[]) -> int:
        """
        Return a random intergene specific breakpoint coordinate.

        Parameters
        ----------
        exclude : List[int]
            list of indices into `self.intergenes`, indicating intergenes that
            should be excluded
            (assume that the coordinates are listed from left to right and
             are contiguous)
        """
        if exclude:
            if exclude[0] < exclude[-1]:        # no wrap
                assert exclude[-1] - exclude[0] == len(exclude) - 1
                range1 = range2 = (0, 0)
                if exclude[0]:
                    range1 = (0, self.intergenes[exclude[0]].sc1 - 1)

                if exclude[-1] != len(self.intergenes) - 1:
                    range2 = (self.intergenes[exclude[-1]].sc2 + 1,
                              self.intergenes[-1].sc2)


                if sum([range1[1], range2[1] - range2[0]]) == 0:
                    raise(CoordinateChoiceError)

                r = G_RNG().choices([range1, range2],
                                  weights=[range1[1], range2[1] - range2[0]])

                return G_RNG().randint(*r[0])
            else:                               # wrap around

                try:
                    return G_RNG().randint(self.intergenes[exclude[-1]].sc2 + 1,
                                           self.intergenes[exclude[0]].sc1 - 1)
                except:
                    raise(CoordinateChoiceError)

        t = sum([x.length for x in self.intergenes]) + len(self.intergenes) - 1
        return G_RNG().randint(0, t)


    def select_random_intergenic_coordinate_excluding(self, c1: int, c2: int,
                                                      d: T_DIR) -> Union[int, None]:
        """
        Return a random intergenic specific coordinate that does not come from
        the intergenes inside of the interval specified by `c1`, `c2`, and `d`.
        The intergenes that include `c1` and `c2` will not contain the returned
        coordinate.

        Parameters
        ----------
        c1 : int
            first intergene specific coordinate
        c2 : int
            second intergene specific coordinate
        d : T_DIR
            the direction {LEFT, RIGHT} to go from `c1`
        """
        try:
            r = self.get_affected_region(c1, c2, d)

            igpositions = r[1]
            c3 = self.select_random_coordinate_in_intergenic_regions(igpositions)

        except CoordinateChoiceError:
            return None

        return c3

    def return_total_coordinate_from_specific_coordinate(self, c, type = INV, debug = False) -> int:

        tc = None
        for r in self.map_of_locations:
            if debug == True:
                print(r)
            if r.itype != type:
                continue
            if r.inSpecific(c):
                distance_to_lower_bound = c - r.sc1
                tc = r.tc1 + distance_to_lower_bound

        if tc is None:
            raise(Exception('Given specific coordinate {c} not found!'))

        return tc

    def return_specific_coordinate_from_total_coordinate(self, c, debug = False):

        sc = None
        for r in self.map_of_locations:
            if debug:
                print(r)
            if r.itype == INV and r.inTotal(c):

                distance_to_lower_bound = c - r.tc1
                sc = r.sc1 + distance_to_lower_bound

                if debug == True:
                    print("PRINTING R")
                    print(r)

        if sc is None:
            raise(Exception('Given total coordinate {c} not found!'))

        return sc


    def get_location_from_coord(self, c: int, use_intergene_specific = False) \
        -> tuple[Interval, int]:
        """
        Given a coordinate, return the endpoints of the Gene or Intergene that
        contains it. Whether total or specific coordinates are used depends
        on the value of `use_intergene_specific`.

        Parameters
        ----------
        c : int
            the coordinate
        use_intergene_specific : bool, optional
            search intergene specific coordinates, if True. otherwise, search
            genes and intergenes using total coordinates.

        Returns
        -------
        Interval
            (interval, gene) where `interval` is Location information for the
            given coordinate, and `gene` is the gene-order position of the
            first occurring gene that occurs at coordinate `c` or after (or the
            length of `self.genes` if no gene is found).
        """
        if not use_intergene_specific:
            for i, l in enumerate(self.map_of_locations):
                if l.inTotal(c):
                    coord = self.return_specific_coordinate_from_total_coordinate(c)
                    gene = self.get_next_gene_location(i)

                    return Interval(*l.asTuple(), c, coord, gene=l.gene), gene
        else:
            for i, l in enumerate(self.map_of_locations):
                if l.isIntergenic() and l.inSpecific(c):
                    coord = self.return_total_coordinate_from_specific_coordinate(c)
                    gene = self.get_next_gene_location(i)

                    return Interval(*l.asTuple(), coord, c, gene=l.gene), gene

        raise(Exception(f'no gene or intergene found for coordinate {c}!'))


    def get_next_gene_location(self, index: int) -> int:
        """
        Given an index into `self.map_of_locations`, return the position of
        the first occurring gene starting from that index. Return the length
        of `self.genes` if no gene is found.

        Parameters
        ----------
        index : int
            index into `self.map_of_locations`
        """
        for j in range(index, len(self.map_of_locations)):
            if self.map_of_locations[j].isGene():
                break

        if not self.map_of_locations[j].isGene():
            return len(self.genes)

        return self.genes.index(cast(Gene, self.map_of_locations[j].gene))


    def return_intergene_by_coordinate(self, c: int) -> Intergene:
        """
        Given a specific coordinate, return the Intergene containing that
        coordinate.
        """

        for intergene in self.iter_intergenes():
            cut1, cut2 = intergene.specific_flanking
            if c >= cut1 and c <= cut2:
                return intergene
        #print(f"Last coordinate is {cut2}")
        raise(Exception(f'no gene or intergene found for coordinate {c}!'))


    def get_affected_region(self, c1: int, c2: int, direction: T_DIR) \
        -> tuple[list[int], list[int], T_PAIR, T_PAIR, Interval, Interval]:
        """
        Return information about the genes and intergenes between the given
        coordinates.

        Parameters
        ----------
        c1 : int
            first intergene specific coordinate
        c2 : int
            second intergene specific coordinate
        direction : T_DIR
            one of LEFT or RIGHT defining whether the region goes left or
            right from c1

        Returns
        -------
        Tuple[List[int], List[int], T_PAIR, T_PAIR, Interval, Interval]
            Returns a tuple
            (genepositions, intergenepositions, firstlengths, secondlengths,
             interval1, interval2)
            where

            0. List of the position of the genes affected. ALWAYS FROM LEFT TO
               RIGHT no matter the `direction`.

            1. List of the position of the intergenes affected. ALWAYS FROM LEFT
               TO RIGHT no matter the `direction`.
               (The intergenes containing c1 and c2 are affected)

            2. Pair with left and right lengths (in nucs) of c1 intergene.
               Watch out, the fact of calling it left or right can be confusing!

            3. Pair with left and right lengths (in nucs) of c2 intergene.
               Same note as above

            4. The intergenic interval containing c1

            5. The intergenic interval containing c2
        """
        l1, _ = self.get_location_from_coord(c1, use_intergene_specific=True)
        l2, _ = self.get_location_from_coord(c2, use_intergene_specific=True)

        p1 = l1.position
        p2 = l2.position

        affected_genes = list()
        affected_intergenes = list()

        t_length = len(self.intergenes)

        if c1 == c2:
            raise(CoordinateChoiceError)

        elif p1 == p2:
            raise(CoordinateChoiceError)

        elif c1 < c2 and direction == T_DIR.RIGHT:

            affected_genes = [i + 1 for i in range(p1, p2)]
            affected_intergenes = [i for i in range(p1, p2 + 1)]

        elif c1 > c2 and direction == T_DIR.RIGHT:

            affected_genes = [i + 1 for i in range(p1, t_length - 1)]
            affected_genes += [i for i in range(0, p2 + 1)]

            affected_intergenes = [i for i in range(p1, t_length)]
            affected_intergenes += [i for i in range(0, p2 + 1)]

        elif c1 > c2 and direction == T_DIR.LEFT:

            affected_genes = [i for i in range(p1, p2, - 1)]
            affected_intergenes = [i for i in range(p1, p2 - 1, -1)]

            affected_genes.reverse()
            affected_intergenes.reverse()

        elif c1 < c2 and direction == T_DIR.LEFT:

            affected_genes = [i for i in range(p1, -1, - 1)]
            affected_genes += [i for i in range(t_length - 1, p2, -1)]

            affected_intergenes = [i for i in range(p1, -1, -1)]
            affected_intergenes += [i for i in range(t_length - 1, p2 - 1, -1)]

            affected_genes.reverse()
            affected_intergenes.reverse()

        c1_lengths = (c1 - l1.sc1, l1.sc2 - c1)
        c2_lengths = (c2 - l2.sc1, l2.sc2 - c2)

        return (affected_genes, affected_intergenes,
                c1_lengths, c2_lengths, l1, l2)
                

    def iter_intergenes(self):
        for x in self.intergenes:
            yield x

    def select_random_length(self, p):
        try:
            return int(af.obtain_value(p, G_NPRNG()))
        except ValueError as e:
            raise ValueError(f"selecting random length with parameter {p}, "
                             f"make sure the number is between 0 and 1: {e}")

    def return_rates(self):

        raise(NotImplementedError)
        # NOT WORKING FOR NOW

        self.total_rates = 0.0

        for gene in self.genes:
            d = gene.family.rates["DUPLICATION"] ## GENE.GENE_FAMILY SHOULD POINT TO A GENE FAMILY OBJECT, NOT A STR
            t = gene.family.rates["TRANSFER"]
            l = gene.family.rates["LOSS"]

            self.total_rates += d + t +l

        return self.total_rates


    def get_num_nucleotides(self) -> int:
        """
        Return the total number of nuclotides in this genome.

        Returns
        -------
        int
            number of nucleotides
        """
        if self._num_nucleotides:
            return self._num_nucleotides
        else:
            assert self.map_of_locations, "map_of_locations not set yet!"
            return self.map_of_locations[-1].tc2


    def get_num_genes(self) -> int:
        """ Return the number of genes in this chromosome. """
        return len(self.genes)

    def get_num_intergenes(self) -> int:
        """ Return the number of intergenes in this chromosome. """
        return len(self.intergenes)

    def __len__(self):
        """ Length of the chromosome in genes. """
        return len(self.genes)

    def __str__(self):

        if self.has_intergenes == True:

            #return ";".join(["CHROMOSOME"] + [str(self.genes[i])+";"+str(self.intergenes[i]) for i in range(len(self.genes))])
            return ";".join([str(self.genes[i]) + ";" + str(self.intergenes[i]) for i in range(len(self.genes))])

        else:

            return ";".join([f"CHROMOSOME {self.name}"] + [str(gene) for gene in self.genes])

    def __iter__(self):

        for x in self.genes:
            yield x

    @abc.abstractmethod
    def obtain_segment(self, gpositions) -> list[Gene]:
        raise(NotImplementedError)

    @abc.abstractmethod
    def invert_segment(self, gpositions):
        raise(NotImplementedError)

    @abc.abstractmethod
    def insert_segment(self, position, segment):
        raise(NotImplementedError)

    @abc.abstractmethod
    def replace_segment(self, affected_positions: list[int],
                        new_segment: list[Gene]):
        raise(NotImplementedError)

    @abc.abstractmethod
    def obtain_intergenic_segment(self, affected_intergenes):
        raise(NotImplementedError)

    @abc.abstractmethod
    def remove_segment(self, segment):
        raise(NotImplementedError)

    @abc.abstractmethod
    def obtain_affected_indices(self, p_extension) -> list[int]:
        raise(NotImplementedError)

    @abc.abstractmethod
    def obtain_affected_indices_family_rates(self, p_extension, gene_families,
                                             mrate) -> list[int]:
        raise(NotImplementedError)

    @abc.abstractmethod
    def cut_and_paste(self, segment: list[Gene],
                      affected_indices: list[int],
                      atposition: int|None=None) -> tuple[int, int]:
        raise(NotImplementedError)

    @abc.abstractmethod
    def obtain_affected_genes_accounting_for_connectedness(self, p_extension,
                                                           interactome):
        raise(NotImplementedError)

    @abc.abstractmethod
    def get_homologous_position(self, segment) -> list[tuple[str, Iterable[int]]]:
        raise(NotImplementedError)


#ChromosomeType = TypeVar('ChromosomeType', bound=Chromosome)
class CircularChromosome(Chromosome):

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            name of the chromosome
        num_nucleotides : int
            length of the chromosome in nucleotides, optional
        """
        super().__init__(*args, **kwargs)
        self.shape = "C"

    def obtain_segment(self, affected_indices: list[int]) -> list[Gene]:
        """
        Get the genes affected by the segment. Return them in the order in
        which they appear in the chromosome (if the segment wraps around, start
        with the rightmost segment).

        Notes
        -----
        We assume that `affected_indices` are given in the order in which
        they appear in the chromosome (i.e. from left to right). So segment
        that wraps will have indices like [8, 9, 0, 1, 2].
        """
        return [self.genes[x] for x in affected_indices]
        #THIS CODE IS UNNECESSARY BECAUSE OF THE ABOVE NOTE:
        #if ends := af.affected_indices_wrap(affected_indices):
        #    return ([self.genes[x] for x in range(ends[0], len(self.genes))] +
        #            [self.genes[x] for x in range(0, ends[1] + 1)])
        #else:
        #    return [self.genes[x] for x in affected_indices]


    def obtain_intergenic_segment(self, affected_intergenes: list[int]) -> list[Intergene]:
        """
        Get the intergenes affected by the segment. Return them in the order in
        which they appear in the chromosome (if the segment wraps around, start
        with the rightmost segment).

        Notes
        -----
        We assume that `affected_indices` are given in the order in which
        they appear in the chromosome (i.e. from left to right). So segment
        that wraps will have indices like [8, 9, 0, 1, 2].
        """
        return [self.intergenes[x] for x in affected_intergenes]
        #if ends := af.affected_indices_wrap(affected_intergenes):
        #    return ([self.intergenes[x] for x in range(ends[0], len(self.intergenes))] +
        #            [self.intergenes[x] for x in range(0, ends[1] + 1)])
        #else:
        #    return [self.intergenes[x] for x in affected_intergenes]


    def remove_segment(self, segment: list[Gene]):
        """ Remove these genes (works with circular chromosomes) """
        for gene in segment:
            self.genes.remove(gene)

    def replace_segment(self, affected_positions: list[int],
                        new_segment: list[Gene]):
        """ Replace the genes in `self.genes` with the new segment.  """
        assert len(affected_positions) == len(new_segment)
        for i, j in enumerate(affected_positions):
            self.genes[j] = new_segment[i]

    def remove_intersegment(self, intersegment: list[Intergene]):

        for intergene in intersegment:
            self.intergenes.remove(intergene)

    def insert_segment(self, position: int, segment: list[Gene]):
        """ Splice the segment into the genes list at this position. """
        self.genes[position:position] = segment


    def invert_segment(self, affected_indices: list[int],
                       affected_inter_indices: list[int]=[]):
        """
        Invert the genes in `self.genes`. Invert all but the first and last
        intergenes in `self.intergenes` if `affected_intergenes` is provided.

        Parameters
        ----------
        affected_genes : List[int]
            the indices of genes to be inverted
        affected_intergenes : List[int]
            the indices of intergenes to be inverted
        """
        reversed_segment = [self.genes[x] for x in reversed(affected_indices)]

        for gene in reversed_segment:
            gene.change_orientation()

        for i, x in enumerate(affected_indices):
            self.genes[x] = reversed_segment[i]

        if affected_inter_indices:         #Now reverse the intergenes:
            intergenes_reversed = [self.intergenes[x]
                                   for x in reversed(affected_inter_indices[1:-1])]

            for revi, i in enumerate(affected_inter_indices[1:-1]):
                self.intergenes[i] = intergenes_reversed[revi]
    

    def invert_divisions(self, cut1, cut2):
        """
        Invert the divisions between the two cuts
        Change the orientation of all divisions affected by the inversion
        """
        start = False
        end = False

        divisions_left = list()
        divisions_right = list()
        divisions_to_invert = list()
        affected_intergenes = list()

        cycle = 0
        
        for intergene in itertools.cycle(self.iter_intergenes()):
            cycle += 1

            if len(intergene.divisions) == 0:
                
                bp = intergene.specific_flanking[0] # breakpoint
                if bp == cut1:
                    start = True
                
                if bp == cut2:
                    end = True

            for division in intergene:

                sf1, sf2 = division.specific_flanking # We get the flankings of the divisions
            
                if cut1 == sf1 and cut2 == sf2:
                    divisions_to_invert.append(division)
                    if intergene not in affected_intergenes:
                        affected_intergenes.append(intergene)
                    end = True
                    break
                
                elif cut1 == sf1 and cut2 != sf2:
                    
                    start = True 
                    divisions_to_invert.append(division)
                    if intergene not in affected_intergenes:
                        affected_intergenes.append(intergene)

                elif cut1 == sf2 and intergene.divisions[-1] == division: # The breakpoint is right before the gene
                    
                    if intergene not in affected_intergenes: # shoudl be ok without this check
                        affected_intergenes.append(intergene)
                    start = True

                elif cut2 == sf1 and start == True:
                    
                    end = True
                    break

                elif cut2 == sf2 and start == True:
                    
                    divisions_to_invert.append(division)
                    if intergene not in affected_intergenes:
                        affected_intergenes.append(intergene)
                    end = True
                    break

                elif start == True:
                    
                    divisions_to_invert.append(division)
                    if intergene not in affected_intergenes:
                        affected_intergenes.append(intergene)
                
            if end == True:
                break    


        all_divisions = list()
    
        for intergene in affected_intergenes:
            for division in intergene:
                all_divisions.append((intergene, division))
       
        #  We prepare the new list with the right order
        
        divisions_left = [division for division in affected_intergenes[0] if division not in divisions_to_invert]
        divisions_right = [division for division in affected_intergenes[-1] if division not in divisions_to_invert]

        right_order_divisions = divisions_left + list(reversed(divisions_to_invert)) + divisions_right

        # We also remove all the divisions in the intergenes affected

        for intergene, division in all_divisions:
            (intergene.divisions).remove(division)

        for division in right_order_divisions:
            for intergene in affected_intergenes:
                empty_space = len(intergene) - sum([len(division) for division in intergene])
                if empty_space >= len(division):
                    if division in divisions_to_invert:
                        division.change_sense()
                    (intergene.divisions).append(division)    
                    break
                    
            

    def inversion_wrap_lengths(self, affected_genes: list[int]
                               ) -> tuple[int, int, int, int]:
        """
        An inversion breaks two breakpoints (in intergenes B1 and B2), replacing
        the gene and intergenes in-place so that, when and inversion wraps
        around, the number of genes and intergenes at the end and beginning are
        not modified.  For example, if we have

        G3 I3 B2 ... B1 G1 I1 G2 I2

        then the result is

        G1 I1 B2 ... B1 G3 I3 G2 I2

        We return the lengths (in breakpoints) of the segments G3 I3 G2 I2 and
        G1 I1. These are the lengths of the segments (at the first part and
        second part) that were inverted, not counting the intergene breakpoints
        in the breakpoint intergenes.

        Notes
        -----
        Assumes the first intergene `self.intergenes[0]` is after the first
        gene `self.genes[0]`.

        Parameters
        ----------
        affected_genes : List[int]
            the indices of the genes affected by an inversion

        Returns
        -------
        Tuple[int, int, int, int]
            (sfirstlen, ssecondlen, tfirstlen, tsecondlen) where sfirstlen
            is the intergene-specific length (in breakpoints) of the right-most
            sequence (not including the B1 intergene), tfirstlen is the total
            length for the same segment. The "second" variables are the same
            values for left-most sequence (not including B2).
        """
        sfirstlen = 0
        tfirstlen = 0
        ssecondlen = 0
        tsecondlen = 0

            #Visit genes and intergenes as pairs 
        igenei = -1
        infirst = True
        for i in range(len(affected_genes)):
            irev = len(affected_genes) - i - 1 
            genei = affected_genes[irev]
            igenei = affected_genes[irev-1]
            if affected_genes[i] == 0:
                infirst = False     # We've wrapped to the beginning
            
            if infirst:
                sfirstlen += self.intergenes[igenei].length + 1
                tfirstlen += self.genes[genei].length + self.intergenes[igenei].length + 1
            else:
                tsecondlen += self.genes[genei].length - 1
                if irev != 0:       # Skip last intergene (it has the breakpoint)
                    ssecondlen += self.intergenes[igenei].length + 1
                    tsecondlen += self.intergenes[igenei].length + 1

        tsecondlen += 1     # Gene breakpoint at index 0 is counted despite
                            # being the same as the intergene breakpoint at -1.
        assert sfirstlen <= tfirstlen and ssecondlen <= tsecondlen
        return sfirstlen, ssecondlen, tfirstlen, tsecondlen


    def cut_and_paste(self, segment: list[Gene],
                      affected_indices: list[int],
                      atposition: int|None=None) -> tuple[int, int]:
        """
        Copy the genes at the given indices to a new location chosen uniformly
        at random.

        Parameters
        ----------
        segment : List[Gene]
            the genes to be moved
        affected_indices : List[int]
            the indices of the genes being moved
        atposition : int, optional
            if provided, the position where to paste the segment (the index
            is in the genome ONCE THE SEGMENT IS REMOVED!), rather than choosing
            a random position.

        Returns
        -------
        tuple[int, int]
            (oldpos, newpos) where `oldpos` is the position to paste the segment
            before the move, and `newpos` is the position in the new chromosome
            where the segment starts
        """
        if len(segment) == len(self.genes):
            return 0, 0

        for gene in segment:
            self.genes.remove(gene)

        if atposition is None:
            position = self.select_random_position()
        else:
            position = atposition

        self.genes[position:position] = segment

        #if ends := af.affected_indices_wrap(affected_indices):
        #    oldpos = position + (ends[1] + 1)
        if affected_indices[0] > affected_indices[-1]:
            oldpos = position + (affected_indices[-1] + 1)
        elif position <= min(affected_indices):
            oldpos = position
        else:
            oldpos = position + len(affected_indices)

        return oldpos, position


    def obtain_affected_indices(self, p_extension) -> list[int]:
        """
        Returns the index list of the affected genes. This will be a range
        of consecutive integers that can WRAP around, but always starts with
        the first affected gene (e.g. [8, 9, 0, 1, 2] when it wraps)
        """
        position = self.select_random_position()
        length = self.select_random_length(p_extension)
        total_length = len(self.genes)

        if length >= total_length:
            return list(range(total_length))

        affected_genes = list()
        for i in range(position, position + length):
            if i >= total_length:
                affected_genes.append(i - total_length)
            else:
                affected_genes.append(i)

        return affected_genes


    def obtain_affected_indices_family_rates(self, p_extension: float,
                                             gene_families: dict[str, GeneFamily],
                                             mrate: str) -> list[int]:

        # In this first version, length is 1. For a more advanced version,
        # I should extent the interactome model
        # Returns N genes accounting for the family rates


        gene2rate = {gene: gene_families[gene.family].rates[mrate] for gene in self.genes}

        if p_extension == 1:

            norm = af.normalize([vl for x, vl in gene2rate.items()])
            mgenes = [i for i in range(len(self.genes))]
            affected_genes = cast(list[int], G_NPRNG().choice(mgenes, size=1, p=norm))
            return affected_genes

        else:

            # If there is an extension

            length = self.select_random_length(p_extension)
            total_length = len(self.genes)

            if length >= total_length:
                affected_genes = [x for x in range(total_length)]
                return affected_genes

            else:
                all_weights = list()

                # If the extension is shorter than the whole genome length

                for position in range(len(self.genes)):
                    affected_genes = list()

                    for i in range(position, position + length):
                        if i >= total_length:
                            affected_genes.append(i - total_length)
                        else:
                            affected_genes.append(i)

                    product = np.prod([gene2rate[self.genes[x]] for x in affected_genes])
                    if product == 0.0:
                        all_weights.append(sys.float_info.min)
                    else:
                        all_weights.append(product)

                position = G_NPRNG().choice(range(len(self.genes)), 1,
                                            p=af.normalize(all_weights))[0]
                affected_genes = list()

                # Returns the index list of the affected genes

                for i in range(position, position + length):
                    if i >= total_length:
                        affected_genes.append(i - total_length)
                    else:
                        affected_genes.append(i)


                return affected_genes


    def obtain_affected_genes_accounting_for_connectedness(self, p_extension,
                                                           interactome: nx.Graph) -> list[int]:
        # Returns N genes accounting for the inverse of the connectedness

        node_degrees = {n:d for n,d in interactome.degree()}    #type: ignore

        corrected_node_degrees = list()
        all_weights = list()

        for gene in self.genes:

            corrected_node_degrees.append(1/(node_degrees[str(gene)] + 1))

        length = self.select_random_length(p_extension)
        total_length = len(self.genes)

        for each_start, gene in enumerate(self.genes):
            position = each_start
            affected_indices = list()

            if length >= total_length:
                # We select the whole genome
                affected_indices = [x for x in range(total_length)]
                return affected_indices

            else:

                for i in range(position, position + length):
                    if i >= total_length:
                        affected_indices.append(i - total_length)
                    else:
                        affected_indices.append(i)

                # We obtain the total weight of this option
                # This means, multiplying all the weights if we start in a given position

                all_weights.append(reduce(lambda x, y: x * y, [corrected_node_degrees[x]
                                                               for x in affected_indices]))

        position = G_NPRNG().choice(range(len(self.genes)), 1,
                                    p=af.normalize(all_weights))[0]
        affected_indices = list()

        # Returns the index list of the affected genes

        for i in range(position, position + length):
            if i >= total_length:
                affected_indices.append(i - total_length)
            else:
                affected_indices.append(i)
        return affected_indices


    def get_homologous_position(self, segment) -> list[tuple[str, Iterable[int]]]:

        homologous = list()

        segment_length = len(segment)
        genes_length = len(self.genes)

        genes = [x.family + "_" + x.orientation for x in self.genes]
        mysegment = [x.family + "_" + x.orientation for x in segment]


        # First we traverse the genome forwards

        for i, gene in enumerate(genes):

            length_counter = 0
            positions: list[int] = list()

            name_gene_in_genome = gene
            name_gene_in_segment = mysegment[0]

            if name_gene_in_genome == name_gene_in_segment:

                positions.append(i)
                length_counter += 1

                for j, x in enumerate(mysegment):
                    if length_counter == segment_length:
                        homologous.append((FORWARDS, tuple(positions)))
                        break
                    if 1 + i + j >= genes_length:
                        if genes[(i + j + 1) - genes_length] == mysegment[j + 1]:
                            positions.append(i+j+1 - genes_length)
                            length_counter += 1
                        else:
                            break
                    else:
                        if genes[i + j + 1] == mysegment[j + 1]:
                            positions.append(i + j + 1)
                            length_counter += 1
                        else:
                            break

        # Second we traverse the genome backwards

        inverted_segment = list()

        for gene in reversed(mysegment):

            inverted_segment.append(gene.replace("+","A").replace("-", "+").replace("A", "-"))

        for i, gene in enumerate(genes):

            positions = list()
            length_counter = 0
            name_gene_in_segment = inverted_segment[0]

            if genes[i] == name_gene_in_segment:

                positions.append(i)
                length_counter += 1

                for j, x in enumerate(inverted_segment):
                    if length_counter == segment_length:
                        homologous.append((BACKWARDS, tuple(positions)))
                        break
                    if 1 + i + j >= genes_length:
                        if genes[(i + j + 1) - genes_length] == inverted_segment[j + 1]:
                            positions.append(i + j + 1 - genes_length)
                            length_counter += 1
                        else:
                            break
                    else:
                        if genes[i + j + 1] == inverted_segment[j + 1]:
                            positions.append(i + j + 1)
                            length_counter += 1
                        else:
                            break

        return homologous

    ### From here, functions related to intergenic regions

    def remove_segment_with_intergenic(self, segment):

        for gene in segment:
            self.genes.remove(gene)

    def insert_gene_within_intergene(self, coordinate, location: Interval, gene):

        tc1, tc2, sc1, sc2, position, t = location.asTuple()

        # Convert to ints:

        #tc1, tc2, sc1, sc2, position = map(int,(tc1,tc2,sc1,sc2,position))

        # The first part is easier - We simply add the gene to the list of genes

        self.genes.insert(position + 1, gene)

        # The second part is cutting the intergene, obtaining the distances

        left_limits = coordinate - sc1
        right_limits = sc2 - coordinate

        # Now we insert the new intergene in the position i + 1

        intergene = Intergene()
        intergene.length = right_limits
        self.intergenes[position].length = left_limits
        self.intergenes.insert(position + 1, intergene)


class LinearChromosome(Chromosome):

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        name : str
            name of the chromosome
        num_nucleotides : int
            length of the chromosome in nucleotides, optional
        """
        super().__init__(*args, **kwargs)
        self.shape = "L"
        raise(NotImplementedError)
    def obtain_segment(self, gpositions) -> list[Gene]:
        raise(NotImplementedError)
    def invert_segment(self, gpositions):
        raise(NotImplementedError)
    def insert_segment(self, position, segment):
        raise(NotImplementedError)
    def obtain_intergenic_segment(self, affected_intergenes):
        raise(NotImplementedError)
    def remove_segment(self, segment):
        raise(NotImplementedError)
    def obtain_affected_indices(self, p_extension) -> list[int]:
        raise(NotImplementedError)
    def obtain_affected_indices_family_rates(self, p_extension, gene_families,
                                             mrate) -> list[int]:
        raise(NotImplementedError)
    def cut_and_paste(self, segment: list[Gene],
                      affected_indices: list[int],
                      atposition: int|None=None) -> tuple[int, int]:
        raise(NotImplementedError)
    def obtain_affected_genes_accounting_for_connectedness(self, p_extension,
                                                           interactome):
        raise(NotImplementedError)
    def get_homologous_position(self, segment) -> list[tuple[str, Iterable[int]]]:
        raise(NotImplementedError)


class Genome:
    """
    Attributes
    ----------
    species: str
        the string indicating the pendant node name
    chromosomes: List[Chromosome]
        the list of chromosomes in this genome
    time: float
        the time in the simulation where this genome was created
        (the time of the most recent event leading to this genome)
    """

    def __init__(self, time=0.0):

        self.species = ""
        self.chromosomes: list[Chromosome] = []
        self.time = time

    def start_genome(self, input):

        for size, shape in input:
            if shape == "L":
                raise(NotImplementedError("Linear chromosomes not implemented"))
                self.chromosomes.append(LinearChromosome("0", size))
            elif shape == "C":
                self.chromosomes.append(CircularChromosome("0", size))

    def select_random_chromosome(self) -> Chromosome:
        """
        At the moment, this just returns the one and only chromosome.
        """
        # I have to weight by the length of each chromosome

        #chromosome = G_NPPRNG.choice(self.chromosomes, 1, p=af.normalize([len(x) for x in self.chromosomes]))[0]

        # So far, only one chromosome per genome, I can safely return the first chromosome
        return self.chromosomes[0]

    def update_genome_species(self, species):

        self.species = species

        for ch in self.chromosomes:
            for gene in ch:
                gene.species = species

    def create_interactome(self, network_model = "BA"):

        self.interactome: nx.Graph
        if network_model == "BA":
            self.interactome  = nx.barabasi_albert_graph(len(self.chromosomes[0]), 1)
        else:
            self.interactome = nx.barabasi_albert_graph(len(self.chromosomes[0]), 1)

        ## Need to shuffle the nodes!!

        randomly_ordered_genes = list(self.chromosomes[0].genes)
        G_RNG().shuffle(randomly_ordered_genes)

        self.interactome = nx.relabel_nodes(self.interactome, {i:str(n) for i,n in enumerate(randomly_ordered_genes)})


    def update_locations(self):
        """ Update the gene and intergene locations in all chromosomes. """
        for ch in self.chromosomes:
            ch.update_locations()


    def init_pieces(self):
        """ Initialize the pieces for all chromosomes. """
        for ch in self.chromosomes:
            ch.fill_pieces()
            ch.update_coordinates()

    #def iter_geneorder_events(self) -> Iterable[GeneOrderEvent]:
    #    """ Get the events from all chromosomes in this genome. """
    #    for chromosome in self.chromosomes:
    #        for event in chromosome.geneorder_history:
    #            yield event

    def __str__(self):

        return ";".join(["GENOME:"] + [str(x) for x in self.chromosomes])

    def __iter__(self):
        for chromosome in self.chromosomes:
            yield chromosome

    def __len__(self):
        return len(self.chromosomes)


class CoordinateChoiceError(Exception):
    pass
