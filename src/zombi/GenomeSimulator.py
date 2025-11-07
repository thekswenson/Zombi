"""
Simulate gene order evolution.

In the Gf mode, genes and intergenic regions are represented explicitly,
knowing both their the intergene-specific coordinates, along with their
"total" genome coordinates. Histories are computed in two steps:
 1. successive random events are applied to the gene/intergene order, and then
 2. the same event order is reconstructed on the resulting "pieces" (genes and
    "divisions", which are subdivisions of intergenes).
"""
import sys
import numpy
import copy
import os
import networkx as nx
import ete3
import itertools

from ete3.coretype.tree import TreeNode
from Bio.SeqFeature import SeqFeature
from typing import Any, Iterable, Union
from pathlib import Path
from collections import Counter, defaultdict

from . import AuxiliarFunctions as af
from .Events import GeneOrderEvent, Loss, Origination, TandemDup, Inversion
from .Events import Transfer, Transposition, MapPseudogeneError, AFER, LFER
from .Events import TDUP, DUP, FER, LOSS, INV, POS, ORIG, GenomeCoordEvent
from .Events import EventOneBreakpoint, EventThreeBreakpoints, LFER_B, LFER_F
from .Events import FORWARDS, BACKWARDS, EventTwoBreakpoints
from .Genomes import Chromosome, CircularChromosome, CoordinateChoiceError
from .Genomes import Gene, GeneFamily, Genome, DivisionFamily, Intergene
from .Genomes import Division, T_DIR, Intergene, LinearChromosome
from .Random import G_RNG, G_NPRNG
from .Filenames import BRANCHEVENTSTABLE, COMPLETETREE, EVENTRATES, FAMILYRATES, GENEORDEREVENTSSCALEDsuffix
from .Filenames import GENEFAMILYGFF, BRANCHEVENTSSCALEDsuffix
from .Filenames import BRANCHEVENTSsuffix, GENEFAMEVENTSsuffix, GENOMEsuffix
from .Filenames import INTERACTOMEsuffix, TREELENGTHS, GENEORDEREVENTSsuffix
from .Filenames import PIECESsuffix, GENEFAMILYLENGTHS, PROFILES
from .Filenames import DIVISIONLENGTHS, LENGTHSsuffix, INITIALGENOMEINFO
from .Filenames import GENEFAMILYINFO, EXTENSIONRATES, TRANSFERRATES


class GenomeSimulator():
    """
    Simulate the evolution of genomes along a species tree.

    Notes
    -----
    Chromosomes have genes and intergenes. Intergenes are divided into
    divisions.  Genes and divisions are both called pieces.
    
    Attributes
    ----------
    active_genomes: set[str]
        list of lineages that currently exist in the gene tree, each lineage is
        represented as a string indicating the pendant node (e.g. n37)
    node_genomes: Dict[str, Genome]
        map lineage name to Genome (gene order) at the nodes of the species
        tree.  this is updated as random events are chosen.
    all_genomes: Dict[str, Genome]
        map lineage name and sequence number (e.g. "n7-3" for the 4th genome on 
        the n7 lineage) to Genome (gene order).
    node_genomes_divs: Dict[str, Genome]
        map lineage name to Genome (gene order) at the nodes of the species
        tree.  this maintains the Divisions in the Intergenes.
    all_genomes_divs: Dict[str, Genome]
        the same as `node_genomes_divs` but for every simulated genome.
        genomes are referenced by lineage name and sequence number (e.g. "n7-3")
    initial_divisions: List[Tuple[int, int]]
        list of tuples (start, end) indicating the divisions in the initial
        genome. For debugging purposes.
    geneorder_events: Dict[str, List[GeneOrderEvent]]
        map lineage name to list of gene order events that have modified a
        genome along that lineage.
    """

    def __init__(self, parameters: dict[str, Any], events_file: Path,
                 root_genome: Path|None, save_all=False):
        """
        Initialize the genome simulator.
        
        Parameters
        ----------
        parameters : dict[str, Any]
            dictionary of parameters for the simulation
        events_file : str
            path to the events file we will read
        root_genome : Path|None
            path to the root genome file in GFF format
        saveall : bool, optional
            whether to save all intermediate genomes, or only those at the
            nodes of the species tree
        """
        self.parameters = parameters

        self.tree_events = self._read_events_file(events_file)
        self.distances_to_start = self._read_distances_to_start(events_file) # Only useful when computing assortative transfers
        self.complete_tree = self._read_tree(events_file.parent / COMPLETETREE)

        self.node_genomes: dict[str, Genome] = dict()
        self.all_genomes: dict[str, Genome] = dict()
        self.all_gene_families: dict[str, GeneFamily] = dict()

        self.all_division_families: dict[int, DivisionFamily] = dict()

        self.gene_families_counter = 0
        self.active_genomes: set[str] = set()
        self.geneorder_events: dict[str, list[GeneOrderEvent]] = defaultdict(list)

        try:
            if self.parameters["RATE_FILE"] != "False":
                if self.parameters["SCALE_RATES"]:
                    self.crown_length = self._read_crown_length(events_file.parent / TREELENGTHS)
                    self.empirical_rates = af.read_empirical_rates(rates_file=self.parameters["RATE_FILE"], scale_rates=self.crown_length)
                else:
                    self.empirical_rates = af.read_empirical_rates(rates_file=self.parameters["RATE_FILE"])
        except KeyError as e:
            sys.exit(f"ERROR: missing parameter {e}.\n"
                     f"       Did you use the correct parameter file?")

        self.root_genome_file = root_genome     #Get root genome from GFF file.
        if root_genome and not root_genome.exists():
            raise(Exception(f"Root genome file {root_genome} not found."))

        # A list to keep track of all the event coordinates
        self.event_coordinates = list() #NOTE: remove this

        # Saving all of the genomes, not just those for the species tree nodes
        self.save_all = save_all
        self.lineage_counter = Counter() #:Current number of saved genomes for this lineage


    def write_genomes(self, genome_folder: Path, intergenic_sequences=False,
                      all_genomes=False):
        """
        Write the genomes to TSV files.
        
        Parameters
        ----------
        intergenic_sequences : bool, optional
            whether to write intergenic sequences lengths, by default False
        allgenomes : bool, optional
            whether to write all simulated genomes, or only those at the
            species tree nodes, by default False
        """
        genome_folder.mkdir(parents=True, exist_ok=True)

        if all_genomes:
            assert self.save_all, "All genomes not saved during the simulation!"
            genomes_to_write = self.all_genomes

        else:
            genomes_to_write = self.node_genomes

        for genome_name, genome in genomes_to_write.items():
            with open(genome_folder / (genome_name + GENOMEsuffix), "w") as f:
                header = ["POSITION", "GENE_FAMILY", "ORIENTATION", "GENE_ID"]
                header = "\t".join(map(str, header)) + "\n"
                f.write(header)

                for chromosome in genome:
                    for index, gene in enumerate(chromosome):

                        line = [index, gene.family, gene.orientation, gene.gene_id]
                        line = "\t".join(map(str,line)) +"\n"
                        f.write(line)

            if intergenic_sequences:
                with open(os.path.join(genome_folder, genome_name + LENGTHSsuffix), "w") as f:
                    header = ["POSITION", "IDENTITY", "LENGTH"]
                    header = "\t".join(map(str, header)) + "\n"
                    f.write(header)

                    for chromosome in genome:
                        i = 0
                        for j, gene in enumerate(chromosome.genes):
                            line = [i, "G(" + str(gene.family) + "_" + str(gene.gene_id) + ")", str(gene.length)]
                            line = "\t".join(map(str, line)) + "\n"
                            f.write(line)
                            i += 1
                            line = [i, INV, str(chromosome.intergenes[j].length)]
                            line = "\t".join(map(str, line)) + "\n"
                            f.write(line)
                            i += 1


    def write_pieces_coordinates(self, genome_folder: Path,
                                 all_genomes: dict[str, Genome]={}):
        """
        Write a TSV file with containing all the pieces of the genome 
        (intergene divisions and genes)
        Genomes have genes and intergenes. Intergenes are divided
        into divisions. Genes and divisions are both called pieces.

        Parameters
        ----------
        all_genomes : bool, optional
            if nonempty, write these genomes instead of self.node_genomes_pieces
        """
        genome_folder.mkdir(parents=True, exist_ok=True)

        if all_genomes:
            genomes_to_write = all_genomes
        else:
            genomes_to_write = self.node_genomes_pieces

        for genome_name, genome in genomes_to_write.items():

            with open(genome_folder / (genome_name + PIECESsuffix), "w") as f:

                header = ["FAMILY", "TYPE", "IDENTITY", "LENGTH", "TOTAL_LEFT",
                          "TOTAL_RIGHT", "ORIENTATION"]
                header = "\t".join(map(str, header)) + "\n"
                f.write(header)

                for chromosome in genome:
                    for piece in chromosome.pieces:
                        if isinstance(piece, Gene):
                            line = "\t".join(list(map(str, [piece.family, piece.ptype, piece.gene_id, piece.length, piece.total_flanking[0], piece.total_flanking[1], piece.orientation ]))) + "\n"
                        else:
                            assert isinstance(piece, Division)
                            line = "\t".join(list(map(str, [piece.family, piece.ptype, piece.identity, piece.length, piece.total_flanking[0], piece.total_flanking[1], piece.orientation ]))) + "\n"
                        f.write(line)


    def write_genome_info(self, genome_folder:Path,
                          filename=INITIALGENOMEINFO):
        """
        Write a TSV file containing gene id, gff gene id, and start and
        end coordinates for every gene and division in the intial genome.
        The coordinates are 1 indexed and inclusive.

        Parameters
        ----------
        genome_folder : str
            the folder
        filename : str, optional
            the filename to use, by default "InitialGenome_info.tsv"
        """
        with open(genome_folder / filename, "w") as f:
            header = ["TYPE", "ID", "START", "END"]
            f.write("\t".join(map(str, header)) + "\n")

            for chromosome in self.initial_genome:
                for gene, intergene in zip(chromosome, chromosome.iter_intergenes()):             
                    line = "\t".join(["GENE_FAMILY", gene.family,
                                      str(gene.start + 1), str(gene.end)]) +\
                                     "\n"
                    f.write(line)
                    end = gene.end
                    for division in intergene:
                        start = end + 1
                        end = start + len(division) - 1
                        line = "\t".join(["DIVISION",
                                          str(division.family),
                                          str(start), str(end)]) + "\n"
                        f.write(line)

    
    def write_gene_family_info(self, genome_folder:Path,
                               filename=GENEFAMILYINFO):
        """
        Write a TSV file with containing gene id, gff gene id, and start and
        end coordinates for every gene in `self.all_gene_families`.

        Parameters
        ----------
        genome_folder : str
            the folder
        filename : str, optional
            the filename to use, by default GENEFAMILYINFO
        """
        with open(genome_folder / filename, "w") as f:
            header = ["GENE_FAMILY", "GFF_ID", "START", "END"]
            f.write("\t".join(map(str, header)) + "\n")

            for gene_family_name, gene_family in self.all_gene_families.items():
                if gene_family.gff_id:
                    line = "\t".join([gene_family_name, str(gene_family.gff_id),
                                      str(gene_family.genes[0].start+1),
                                      str(gene_family.genes[0].end)]) + "\n"
                    f.write(line)


    def write_gene_family_GFF_ids(self, genome_folder:str,
                               filename=GENEFAMILYGFF):
        """
        Write a TSV file with containing gene id and gff gene id for every gene
        in `self.all_gene_families`.

        Parameters
        ----------
        genome_folder : str
            the folder
        filename : str, optional
            the filename to use, by default GENEFAMILYGFF
        """
        with open(os.path.join(genome_folder, filename), "w") as f:
            header = ["GENE_FAMILY", "GFF_ID"]
            f.write("\t".join(map(str, header)) + "\n")

            for gene_family_name, gene_family in self.all_gene_families.items():
                if gene_family.gff_id:
                    line = "\t".join([gene_family_name, str(gene_family.gff_id)]) + "\n"
                    f.write(line)


    def write_gene_family_lengths(self, genome_folder):

        with open(os.path.join(genome_folder, GENEFAMILYLENGTHS), "w") as f:
            header = ["GENE_FAMILY", "LENGTH"]
            header = "\t".join(map(str, header)) + "\n"
            f.write(header)

            for gene_family_name, gene_family in self.all_gene_families.items():
                line = "\t".join([gene_family_name, str(gene_family.length)]) + "\n"
                f.write(line)


    def write_division_lengths(self, genome_folder):

        with open(os.path.join(genome_folder, DIVISIONLENGTHS), "w") as f:
            header = ["DIVISION_ID", "LENGTH"]
            header = "\t".join(map(str, header)) + "\n"
            f.write(header)

            for div_family_name, div_family in self.all_division_families.items():
                line = "\t".join([str(div_family_name), str(len(div_family))]) + "\n"
                f.write(line)


    def write_gene_family_events(self, gene_family_events_folder: Path):

        if not os.path.isdir(gene_family_events_folder):
            os.mkdir(gene_family_events_folder)

        for gene_family_name, gene_family in self.all_gene_families.items():

            with open(gene_family_events_folder / (gene_family_name + GENEFAMEVENTSsuffix),"w") as f:

                header = ["TIME","EVENT","NODES"]

                header = "\t".join(map(str, header)) + "\n"

                f.write(header)

                for time, event, nodes in gene_family.events:
                    
                    line = [time, event, nodes]

                    line = "\t".join(map(str, line)) + "\n"
                    f.write(line)
        

    def write_gene_trees(self, gene_tree_folder: Path, gene_trees = True, reconciliations = False):

        gene_tree_folder.mkdir(parents=True, exist_ok=True)

        for gene_family_name, gene_family in self.all_gene_families.items():

            complete_tree, pruned_tree, rec_tree = gene_family.generate_tree()

            if gene_trees == True:

                with open(gene_tree_folder / (gene_family_name + "_completetree.nwk"), "w") as f:
                    f.write(complete_tree)
                if pruned_tree != None:
                    with open(gene_tree_folder / (gene_family_name + "_prunedtree.nwk"), "w") as f:
                        f.write(pruned_tree)
            if reconciliations == True:
                with open(gene_tree_folder / (gene_family_name + "_rec.xml"), "w") as f:
                    f.write(rec_tree)


    def write_geneorder_events_per_branch(self, out_dir: Path, scale: float,
                                          scaled_file: Path, events_file: Path):
        """
        Organize the gene order events per branch and write them to files.
        Each branch will have a TSV file with columns:

        TIME    EVENT   BREAKPOINTS  CHROMOSOME

        where BREAKPOINTS are the gene-order positions affected by the event,
        and CHROMOSOME is the name of the chromosome where the event happened.
        """
        out_dir.mkdir(parents=True, exist_ok=True)

        for lineage, events in self.geneorder_events.items():
            write_geneorder_events(out_dir / (lineage + GENEORDEREVENTSsuffix), events)
        if scale:
            for lineage, events in self.geneorder_events.items():
                write_geneorder_events(out_dir / (lineage + GENEORDEREVENTSSCALEDsuffix),
                                       events, scale, scaled_file, events_file)


    def write_events_per_branch(self, events_per_branch_folder: Path,
                                scale: float, scaled_file: Path,
                                events_file: Path):
        ### THIS FUNCTION SHOULD BE CLEANED! THE INFO NOW IS REDUNDANT
        
        def clever_writing():
            table = list()
            for genome_name, genome in self.node_genomes.items():            
                for chromosome in genome:                
                    
                    for event in chromosome.event_history:                                  
                        
                        etype, time, breakpoints = event.return_info()
                        
                        if etype == FER:
                            assert isinstance(event, Transfer)
                            breakpoints = str(event.receptortbp)
                            
                        if etype == POS:
                            assert isinstance(event, Transposition)
                            breakpoints += "," + str(event.tbpH)  
                        
                        table.append((genome_name, time, etype, breakpoints))
                        
            table = sorted(table, key=lambda x:x[1])
            
            with open(events_per_branch_folder / BRANCHEVENTSTABLE, "w") as f:
                
                header = "\t".join(["BRANCH", "TIME", "EVENT", "BREAKPOINTS"]) + "\n"
                f.write(header)
                
                for e in table:
                    f.write("\t".join([str(x) for x in e]) + "\n")
        
        events_per_branch_folder.mkdir(parents=True, exist_ok=True)

        clever_writing()

        # Build events per branch
        events_per_branch = dict()

        for gene_family_name, gene_family in self.all_gene_families.items():

            for time, event, nodes in gene_family.events:
                
                name = nodes.split(";")[0]

                if name not in events_per_branch:
                    events_per_branch[name] = list()

                if event == "S" or event == "E" or event == "F":
                    continue

                elif event == FER:
                    donor = name
                    recipient = nodes.split(";")[4]

                    handle = nodes.split(";")
                    gene_names = list(map(lambda x: gene_family_name + "_" + x, [handle[1], handle[3], handle[5]]))
                    new_nodes = ";".join([donor, gene_names[0], donor, gene_names[1], recipient, gene_names[2]])


                    if donor not in events_per_branch:
                        events_per_branch[donor] = list()

                    events_per_branch[donor].append((time, "LT", new_nodes))

                    if recipient not in events_per_branch:
                        events_per_branch[recipient] = list()

                    events_per_branch[recipient].append((time, "AT", new_nodes))

                elif event == TDUP:

                    handle = nodes.split(";")
                    new_nodes = ";".join(map(lambda x: gene_family_name + "_" + x,[handle[1],handle[3],handle[5]]))
                    events_per_branch[name].append((time, event, new_nodes))

                else:

                    gene_id = nodes.split(";")[-1]
                    events_per_branch[name].append((time, event, gene_family_name + "_" + gene_id))


        for name, events in events_per_branch.items():

            with open(events_per_branch_folder / (name + BRANCHEVENTSsuffix), "w") as f:

                header = ["TIME", "EVENT", "NODES"]
                header = "\t".join(map(str, header)) + "\n"

                f.write(header)

                for time, event, nodes in sorted(events, key = lambda x: float(x[0])):

                    line = [str(time), event, nodes]
                    line = "\t".join(line) + "\n"
                    f.write(line)
        
            if scale: # Only working if Species Tree has been scaled to the same distance!

                beginning_time, mfactor = get_scaling_factor(scaled_file, events_file, scale)

                with open(events_per_branch_folder / (name + BRANCHEVENTSSCALEDsuffix), "w") as f:

                    header = ["TIME", "EVENT", "NODES"]
                    header = "\t".join(map(str, header)) + "\n"

                    f.write(header)

                    for time, event, nodes in sorted(events, key = lambda x: float(x[0])):
                        time = (float(time) - beginning_time) * mfactor
                        line = [str(time), event, nodes]
                        line = "\t".join(line) + "\n"
                        f.write(line)

    
    
    def write_profiles(self, profiles_folder):

        if not os.path.isdir(profiles_folder):
            os.mkdir(profiles_folder)


        genome_names = [x for x in self.node_genomes.keys()]
        gene_family_names = [str(x) for x in self.all_gene_families.keys()]

        # For clarity, I start with Initial Genome

        genome_names[0], genome_names[1] = genome_names[1], genome_names[0]

        data = list()
        data.append(["GENOME"] + gene_family_names)


        for genome_name in genome_names:
            line = dict()

            genome = self.node_genomes[genome_name]

            for gene_family_name in gene_family_names:
                if gene_family_name not in line:
                    line[gene_family_name] = 0

            for chromosome in genome:
                for index, gene in enumerate(chromosome):
                    line[gene.family] += 1

            data.append([genome_name] + [str(line[fm]) for fm in gene_family_names])

        # We will transpose the data

        mlenx = len(data)
        mleny = len(data[0])

        with open(os.path.join(profiles_folder, PROFILES), "w") as f:
            for i in range(mleny):
                line = list()
                for j in range(mlenx):
                    line.append(str(data[j][i]))
                line = "\t".join(line) + "\n"
                f.write(line)


    def write_interactomes(self, genome_folder):

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        for genome_name, genome in self.node_genomes.items():

            if not hasattr(genome, "interactome"):
                continue

            with open(os.path.join(genome_folder, genome_name + INTERACTOMEsuffix), "w") as f:

                header = ["GENE_1", "GENE_2"]
                header = "\t".join(map(str, header)) + "\n"
                f.write(header)

                for g1, g2 in genome.interactome.edges:
                    f.write("\t".join([g1, g2]) + "\n")

    def write_family_rates(self, genome_folder):

        with open(os.path.join(genome_folder, FAMILYRATES), "w") as f:
            header = ["GENE_FAMILY", TDUP, FER, LOSS]
            header = "\t".join(map(str, header)) + "\n"
            f.write(header)

            for gene_family_name, gene_family in self.all_gene_families.items():

                d = gene_family.rates["TANDEMDUP"]
                t = gene_family.rates["TRANSFER"]
                l = gene_family.rates["LOSS"]

                f.write("\t".join(map(str,[gene_family_name, d,t,l])) + "\n")

    def _read_events_genome(self, events_file):

        events = list()
        with open(events_file) as f:
            f.readline()
            for line in f:
                handle = line.strip().split("\t")
                events.append(handle)
        return events
    
    def _read_events_file(self, events_file) -> list[tuple[float, str, str]]:
        """
        Return a list of tuples (time, event, nodes) from the Events.tsv file.
        """
        events = list()
        with open(events_file) as f:
            f.readline()
            for line in f:
                time, event, nodes = line.strip().split("\t")
                events.append((float(time), event, nodes))

        return events
    
    def read_genome_events_file(self, events_file):

        events = list()
        with open(events_file) as f:
            for line in f:
                h = line.strip().split("\t")
                events.append((h[0], float(h[1]), h[2], h[3], (int(h[4]), int(h[5]), False)))
        return events

    def write_event_file(self, event, event_file):

        # For debugging purposes

        line = ["G"]

        with open(event_file, "a") as f:
            line.append(str(event.time))
            line.append(str(event.etype))
            line.append(str(event.lineage))
            line.append(str(event.sbpL))
            line.append(str(event.sbpR))
            line.append("RIGHT")
            line = "\t".join(line) + "\n"
            f.write(line)
            
    def _read_distances_to_start(self, events_file):

        # This function could be fusion with the function above

        distances_to_start = dict()

        with open(events_file) as f:
            f.readline()
            for line in f:
                time, _, nodes = line.strip().split("\t")
                n = nodes.split(";")[0]
                distances_to_start[n] = float(time)
        return distances_to_start


    def _read_crown_length(self, length_file):

        with open(length_file) as f:

            cl = float(f.readlines()[-1].strip().split("\t")[-1])

        return cl

    def _read_tree(self, tree_file):

        with open(tree_file) as f:
            t = ete3.Tree(f.readline().strip(), format=1)
        return t

    def return_new_identifiers_for_segment(self, segment: list[Gene]):
        """ For each gene get the next unused identifier for it's family. """
        new_identifiers = list()

        for gene in segment:
            gf = gene.family
            new_id = self.all_gene_families[gf].obtain_new_gene_id()
            new_identifiers.append(new_id)

        return new_identifiers

    def copy_segment(self, segment: list[Gene]) -> list[Gene]:
        """
        Deep copy the genes of the segment, while updating the gene_id.
        """
        newids = self.return_new_identifiers_for_segment(segment)
    
        new_segment = []
        for gene, newid in zip(segment, newids):
            new_gene = copy.deepcopy(gene)
            new_gene.gene_id = newid
            new_segment.append(new_gene)
    
        return new_segment

    def return_new_gene_ids_for_segment_with_divisions(self, segment: list[Gene | Division]) \
        -> list[int]:

        new_identifiers = list()

        for gene in segment:
            if not isinstance(gene, Gene):
                continue
            gf = gene.family
            new_id = self.gene_families_second[gf].obtain_new_gene_id()
            new_identifiers.append(new_id)

        return new_identifiers


    def fill_genome(self, intergenic_sequences = False, family_rates = False,
                    interactome = False):
        """
        Create the initial genome, which is a single circular chromosome with
        genes and intergenic regions in order.
        """
        genome = Genome()
        genome.species = "Root"
        time = 0

        initial_genome_size = self.parameters["INITIAL_GENOME_SIZE"].split(";")
        shape = "C"

        for chromnum, n_genes in enumerate(initial_genome_size):

            if shape == "L":
                raise(NotImplementedError("Linear chromosomes not implemented"))
                chromosome = LinearChromosome(str(chromnum), 0)
            elif shape == "C":
                chromosome = CircularChromosome(str(chromnum))
            else:
                raise(Exception('unexpected chromosome shape'))

            if intergenic_sequences == True:

                chromosome.has_intergenes = True
                mean_length = int(self.parameters["INTERGENE_LENGTH"])
                intergene_lengths = [int(x * mean_length * int(n_genes)) for x in
                                     af.sample_from_dirichlet(int(n_genes), G_NPRNG())]

                for i in range(int(n_genes)):
                    intergenic_sequence = Intergene()
                    intergenic_sequence.length = intergene_lengths[i]
                    chromosome.intergenes.append(intergenic_sequence)

            for i in range(int(n_genes)):

                # We fill the chromosomes and we create also the gene families
                if family_rates == True and self.parameters["RATE_FILE"] == "False":
                    gene, gene_family = self.make_origination(genome.species, time, family_mode=True)

                elif family_rates == True and self.parameters["RATE_FILE"] != "False":
                    gene, gene_family = self.make_origination(genome.species, time, family_mode=True,
                                                              empirical_rates=True)
                else:
                    gene, gene_family = self.make_origination(genome.species, time)

                initial_gene = copy.deepcopy(gene)
                initial_gene.species = "Initial"

                gene_family.genes.append(initial_gene)
                chromosome.genes.append(gene)

                self.all_gene_families[str(self.gene_families_counter)] = gene_family

                if intergenic_sequences == True:

                    gene.length = int(af.obtain_value(self.parameters["GENE_LENGTH"], G_NPRNG()))

            if intergenic_sequences == True:
                chromosome.update_locations()

            genome.chromosomes.append(chromosome)

            if interactome == True:
                genome.create_interactome()

        self.initial_genome = copy.deepcopy(genome)
        self.initial_gene_families =  copy.deepcopy(self.all_gene_families)

        return genome


    def read_genome(self, genome_file: Path, intergenic_sequences = False,
                    family_rates = False, interactome = False):
        """
        Create a genome with genes and intergenic regions specified by the given
        `genome_file` (.gff).

        Parameters
        ----------
        genome_file : Path
            the filename (.gff) that annotates the genes in the genome. A
            single chromosome is assumed.

        Returns
        -------
        Genome
            the newly constructed genome with gene and intergenic sizes having
            the lengths specified in `genome_file`.

        Notes
        -----
            Multiple chromosomes are not supported at the moment.
        """
        genome = Genome()
        genome.species = "Root"
        time = 0

        chrom_len, gene_features = af.parse_GFF(genome_file)

            #Create a chromosome of the appropriate shape:
        shape = "C"
        chromosome: Chromosome
        if shape == "L":
            raise(NotImplementedError("Linear chromosomes not implemented"))
            chromosome = LinearChromosome("0", chrom_len)
        else:
            assert shape == "C"
            chromosome = CircularChromosome("0", chrom_len)

            #Create the genes:
        prev_feature = None     #previous feature used to make a gene
        for feature in gene_features:
            if prev_feature and prev_feature.location.end > feature.location.start:  #type: ignore
                print(f'WARNING: skipping the creation of overlapping gene '
                      f'{feature.id},\n\tas it overlaps with {prev_feature.id}')
            else:
                gene, gene_family = self.make_gene(feature, genome.species,
                                                   time, family_rates,
                                                   self.parameters["RATE_FILE"] != "False")
                chromosome.genes.append(gene)
                self.all_gene_families[str(self.gene_families_counter)] = gene_family
                gene_family = gene.orientation

                prev_feature = feature

            #Create the intergenes:
        if intergenic_sequences:    #first intergene is after the first gene
            chromosome.has_intergenes = True

            if shape == "L":        #before the first gene
                chromosome.intergenes.append(Intergene(chromosome.genes[0].start))

            for gene1, gene2 in itertools.pairwise(chromosome.genes):
                chromosome.intergenes.append(Intergene(gene2.start - gene1.end))

            if shape == "L":        #after the last gene
                intergene = Intergene(chrom_len - chromosome.genes[-1].end)
            elif shape == "C":      #betweeen the last and first genes
                intergene = Intergene((chrom_len - chromosome.genes[-1].end) +
                                      chromosome.genes[0].start)
            else:
                raise(Exception(f"Unrecognized shape string: {shape}"))
            chromosome.intergenes.append(intergene)

                                #NOTE: this is called in run_f as well!
            chromosome.update_locations()

        genome.chromosomes.append(chromosome)

        if interactome:
            genome.create_interactome()

        self.initial_genome = copy.deepcopy(genome)
        self.initial_gene_families =  copy.deepcopy(self.all_gene_families)

        return genome


    def make_gene(self, gene_feature: SeqFeature, species_tree_node: str,
                  time: int, family_mode = False, empirical_rates = False) -> tuple[Gene, GeneFamily]:
        """
        Make a new gene in a new gene family, based on the given `gene_feature`.

        Parameters
        ----------
        gene_feature : SeqFeature
            the Biopython gene feature
        species_tree_node : str
            the species tree node to assign to this genome
        time : int
            the time tick when this genome exists
        family_mode : bool, optional
            the gene family should have its own rates, by default False
        empirical_rates : bool, optional
            use empirical rates rather than new rates, by default False

        Returns
        -------
        Tuple[Gene, GeneFamily]
            [description]
        """
        self.gene_families_counter += 1
        gene_family_id = str(self.gene_families_counter)

        gene = Gene()

        if gene_feature.strand == 1:
            gene.orientation = "+"
        elif gene_feature.strand == -1:
            gene.orientation = "-"
        else:
            raise(Exception(f"Unknown strand for gene:\n{gene_feature}"))

        gene.family = str(self.gene_families_counter)
        gene.species = species_tree_node
        gene.start = gene_feature.location.start    #type: ignore
        gene.end = gene_feature.location.end        #type: ignore
        gene.length = gene.end - gene.start         #type: ignore

        gene_family = GeneFamily(gene_family_id, time)
        gene_family.length = gene.length
        gene_family.genes.append(gene)
        gene_family.gff_id = gene_feature.id
        gene.gene_id = gene_family.obtain_new_gene_id()

        self.all_gene_families[gene_family_id] = gene_family
        self.all_gene_families[gene.family].append_event(time, ORIG, species_tree_node)

        if family_mode and not empirical_rates:
            d, t,l, _, _, _ = self.generate_new_rates()
            gene_family.rates["TANDEMDUP"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        elif family_mode and empirical_rates:
            d, t, l = self.generate_empirical_rates()
            gene_family.rates["TANDEMDUP"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        initial_gene = copy.deepcopy(gene)
        initial_gene.species = "Initial"
        gene_family.genes.append(initial_gene)

        return gene, gene_family


    def run(self):
        """
        Simulate the genomes with no intergenes.
        """
        d = af.obtain_value(self.parameters["TANDEMDUP"], G_NPRNG())
        u = af.obtain_value(self.parameters["DUPLICATION"], G_NPRNG())
        t = af.obtain_value(self.parameters["TRANSFER"], G_NPRNG())
        l = af.obtain_value(self.parameters["LOSS"], G_NPRNG())
        i = af.obtain_value(self.parameters["INVERSION"], G_NPRNG())
        c = af.obtain_value(self.parameters["TRANSPOSITION"], G_NPRNG())

        o = af.obtain_value(self.parameters["ORIGINATION"], G_NPRNG())

        # First we prepare the first genome

        genome = self.fill_genome()

        self.active_genomes.add("Root")
        self.node_genomes["Root"] = genome

        # We add the original genome too

        self.node_genomes["Initial"] = copy.deepcopy(genome)
        if self.save_all:
            self.all_genomes["Initial-0"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)
        # Second, we compute the time to the next event:

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes),
                                                                    [d, u, t, l, i, c, o])

            elapsed_time = float(current_time) - elapsed_time

            if time_to_next_genome_event + current_time >= float(time_of_next_species_tree_event):

                current_species_tree_event +=1
                current_time = time_of_next_species_tree_event

                if event == "S":
                    sp,c1,c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes
                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes
                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)

                    self.node_genomes[c1] = genome_c1
                    self.node_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:
                current_time += time_to_next_genome_event

                if info := self.evolve_genomes(d, u, t, l, i, c, o, current_time):
                    if self.save_all and info:
                        self.add_genome_to_all_genomes(info[1])


    def run_i(self):
        """
        This is interactive mode, which doesn't appear in the documentation.
        """

        # Interactome mode

        d = af.obtain_value(self.parameters["DUPLICATION"], G_NPRNG())
        t = af.obtain_value(self.parameters["TRANSFER"], G_NPRNG())
        l = af.obtain_value(self.parameters["LOSS"], G_NPRNG())
        i = af.obtain_value(self.parameters["INVERSION"], G_NPRNG())
        c = af.obtain_value(self.parameters["TRANSPOSITION"], G_NPRNG())
        o = af.obtain_value(self.parameters["ORIGINATION"], G_NPRNG())
        rm = af.obtain_value(self.parameters["REMOVE"], G_NPRNG())
        rw = af.obtain_value(self.parameters["REWIRE"], G_NPRNG())

        # First we prepare the first genome

        genome = self.fill_genome(interactome = True)

        # We prepare to important dicts in this mode

        self.active_genomes.add(genome.species)
        self.node_genomes["Root"] = genome

        # We add the initial genome too

        self.node_genomes["Initial"] = copy.deepcopy(genome)
        if self.save_all:
            self.all_genomes["Initial-0"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)

        # Second, we compute the time to the next event:

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), [d, t, l, i, c, o, rm, rw])

            elapsed_time = float(current_time) - elapsed_time

            if time_to_next_genome_event + current_time >= float(time_of_next_species_tree_event):

                current_species_tree_event +=1
                current_time = time_of_next_species_tree_event

                if event == "S":
                    sp,c1,c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes
                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)

                    self.node_genomes[c1] = genome_c1
                    self.node_genomes[c2] = genome_c2


                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                if info := self.evolve_genomes_i(d, t, l, i, c, o, rm, rw, current_time):
                    if self.save_all and info:
                        self.add_genome_to_all_genomes(info[1])


    def run_m(self):

        if self.parameters["RATE_FILE"] == "False":
            sys.exit("Error: Family mode (Gm) requires RATE_FILE to be set.")

        # First we prepare the first genome

        genome = self.fill_genome(family_rates=True)

        # We prepare to important dicts in this mode

        self.active_genomes.add(genome.species)
        self.node_genomes["Root"] = genome

        # We add the initial genome too

        self.node_genomes["Initial"] = copy.deepcopy(genome)
        if self.save_all:
            self.all_genomes["Initial-0"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)

        # Second, we compute the time to the next event:

        elapsed_time = 0.0


        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event_family_mode()

            elapsed_time = float(current_time) - elapsed_time

            if time_to_next_genome_event + current_time >= float(time_of_next_species_tree_event):

                current_species_tree_event +=1
                current_time = time_of_next_species_tree_event

                if event == "S":
                    sp,c1,c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes
                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)

                    self.node_genomes[c1] = genome_c1
                    self.node_genomes[c2] = genome_c2


                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                if info := self.evolve_genomes_m(current_time):
                    if self.save_all and info:
                        self.add_genome_to_all_genomes(info[1])


    def run_u(self):
        """
        Simulate genomes with user defined genome rates (using RatesCustomizer).
        """
        genome = self.fill_genome()
        self.active_genomes.add(genome.species)
        self.node_genomes["Root"] = genome

        # We add the original genome too

        self.node_genomes["Initial"] = copy.deepcopy(genome)
        if self.save_all:
            self.all_genomes["Initial-0"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event_advanced_modes()

            elapsed_time = float(current_time) - elapsed_time
            if time_to_next_genome_event + current_time >= float(time_of_next_species_tree_event):
                current_species_tree_event += 1
                current_time = time_of_next_species_tree_event
                if event == "S":
                    sp, c1, c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes
                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)
                    self.node_genomes[c1] = genome_c1
                    self.node_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                if info := self.advanced_evolve_genomes(current_time):
                    if self.save_all and info:
                        self.add_genome_to_all_genomes(info[1])


    def run_f(self):

        d = af.obtain_value(self.parameters["TANDEMDUP"], G_NPRNG())
        u = af.obtain_value(self.parameters["DUPLICATION"], G_NPRNG())
        t = af.obtain_value(self.parameters["TRANSFER"], G_NPRNG())
        l = af.obtain_value(self.parameters["LOSS"], G_NPRNG())
        i = af.obtain_value(self.parameters["INVERSION"], G_NPRNG())
        c = af.obtain_value(self.parameters["TRANSPOSITION"], G_NPRNG())
        o = af.obtain_value(self.parameters["ORIGINATION"], G_NPRNG())

        # First we prepare the root genome

        if self.root_genome_file:
            genome = self.read_genome(self.root_genome_file,
                                      intergenic_sequences=True)
        else:
            genome = self.fill_genome(intergenic_sequences=True)

        ## These two lines are important for this mode (already in read_genome)

        #for chromosome in genome:          #NOTE: already done in read_genome!?
        #    chromosome.obtain_flankings()
        #    chromosome.obtain_locations()

        self.active_genomes.add(genome.species)
        self.node_genomes["Root"] = genome

        # We add the original genome too

        self.node_genomes["Initial"] = copy.deepcopy(genome)
        if self.save_all:
            self.all_genomes["Initial-0"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)

        # Second, we compute the time to the next event:

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), [d, t, l, i, c, o])

            elapsed_time = float(current_time) - elapsed_time

            if time_to_next_genome_event + current_time >= float(time_of_next_species_tree_event):

                current_species_tree_event += 1
                current_time = time_of_next_species_tree_event

                if event == "S":
                    sp, c1, c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes

                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)

                    self.node_genomes[c1] = genome_c1
                    self.node_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:
                current_time += time_to_next_genome_event
                if info := self.advanced_evolve_genomes_f(d, u, t, l, i, c, o, current_time):
                    if self.save_all and info:
                        self.add_genome_to_all_genomes(info[1])


    def run_f_debug(self, injected_events): # Only for debugging purposes

        # First we prepare the root genome

        if self.root_genome_file:
            genome = self.read_genome(self.root_genome_file,
                                      intergenic_sequences=True)
        else:
            genome = self.fill_genome(intergenic_sequences=True)


        self.active_genomes.add(genome.species)
        self.node_genomes["Root"] = genome

        # We add the original genome too

        self.node_genomes["Initial"] = copy.deepcopy(genome)
        if self.save_all:
            self.all_genomes["Initial-0"] = copy.deepcopy(genome)

        all_species_tree_events = [(FER,float(x[0]),x[1],x[2]) for x in self.tree_events]
        
        # Second, we compute the time to the next event:

        all_events = sorted(all_species_tree_events + injected_events, key=lambda x: x[1]) # We put together the two types of events and order by time
        

        for item in all_events:

            if item[0] == FER:
                etype, time, event, nodes = item

                if event == "S":
                    sp, c1, c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes

                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, time)

                    self.node_genomes[c1] = genome_c1
                    self.node_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(time)
                    break

            else:
                # If the event is a Genome level event
                #print(item)

                etype, time, event, nodes, r = item
                
                if event == POS:
                    c1, c2, c3, d = r
                elif event == ORIG:
                    c1, d = r
                elif event == FER:
                    c1, c2, c3, d, lineage_r = r
                else:
                    c1, c2, d = r

 
                if event == TDUP:

                    lineage = nodes
                    ch = self.node_genomes[lineage].chromosomes[0] 
                    self.update_genome_indices(lineage)

                    self.make_tandemdup_intergenic(ch, c1, c2, d, lineage, time)

                elif event == FER:

                    lineage = nodes
                    ch = self.node_genomes[lineage].chromosomes[0] 
                    chreceptor = self.node_genomes[lineage_r].chromosomes[0] 
                    
                    self.update_genome_indices(lineage)
                    self.update_genome_indices(lineage_r)

                    

                    self.make_transfer_intergenic(ch, c1, c2, d, lineage, chreceptor,
                                                  c3, lineage_r, time)
                    
                elif event == LOSS:

                    lineage = nodes
                    pseudo = False
                    if G_NPRNG().uniform(0,1) <= float(self.parameters["PSEUDOGENIZATION"]):
                        pseudo = True

                    ch = self.node_genomes[lineage].chromosomes[0]
                    self.update_genome_indices(lineage)
                    self.make_loss_intergenic(ch, c1, c2, d, lineage, time, pseudo)
                    

                elif event == INV:

                    lineage = nodes
                    self.update_genome_indices(lineage)
                    ch = self.node_genomes[lineage].chromosomes[0] 

                    self.make_inversion_intergenic(ch, c1, c2, d, lineage, time)

                elif event == POS:

                    lineage = nodes
                    self.update_genome_indices(lineage)
                    ch = self.node_genomes[lineage].chromosomes[0] 
                    
                    self.make_transposition_intergenic(ch, c1, c2, d, c3, lineage, time)
    

                elif event == ORIG:

                    lineage = nodes
                    self.update_genome_indices(lineage)
                    ch = self.node_genomes[lineage].chromosomes[0]
                    self.make_origination_intergenic(ch,c1,lineage, time)


    def add_genome_to_all_genomes(self, toadd: Genome|str|tuple[str, str],
                                  checklineage=False):
        """
        Add a genome or a pair of genomes to the `self.all_genomes` dict.

        Parameters
        ----------
        toadd : Genome | str | tuple[str, str]
            If `toadd` is a Genome, add it to the `all_genomes` dict.
            If it is a string, add the current genome from lineage `toadd`.
            If it is a tuple, add both genomes from the specified lineages.
        checklineage : bool, optional
            only add the genome if the lineage counter for that genome is 0
            (in other words, no genome exists for that lineage)
        """
        if isinstance(toadd, str):
            lineage = toadd
            genome = copy.deepcopy(self.node_genomes[lineage])

        elif isinstance(toadd, tuple):
            assert len(toadd) == 2
            self.add_genome_to_all_genomes(toadd[0], checklineage)
            self.add_genome_to_all_genomes(toadd[1], checklineage)
            return

        else:
            lineage = toadd.species
            genome = copy.deepcopy(toadd)

        if not checklineage or self.lineage_counter[lineage] == 0:
            self.all_genomes[f'{lineage}-{self.lineage_counter[lineage]}'] = genome
            self.lineage_counter[lineage] += 1


    def add_pieces(self, lineage: str, all_genomes_pieces: dict[str, Genome],
                   lineage_counter: Counter[str], checklineage=False):
        """
        Add the current genome in `self.node_genomes_pieces` for the given
        `lineage` to the `all_genomes_pieces` dict. If the counter for
        lineage "n3" is at 2, then the key will be "n3-2".
        
        Parameters
        ----------
        checklineage : bool, optional
            only add the genome if the lineage counter for that genome is 0
            (in other words, no genome exists for that lineage)
        """
        if not checklineage or lineage_counter[lineage] == 0:
            genome = copy.deepcopy(self.node_genomes_pieces[lineage])
            all_genomes_pieces[f'{lineage}-{lineage_counter[lineage]}'] = genome
            lineage_counter[lineage] += 1


    def update_genome_indices(self, lineage):
        """
        Update the indices for genes and intergenes. This should be called after
        every rearrangement (NOTE: why is it not called within the rearrangement
        code?).
        """
        self.node_genomes[lineage].update_locations()
    

    def generate_new_rates(self):

        d = af.obtain_value(self.parameters["DUPLICATION"], G_NPRNG())
        t = af.obtain_value(self.parameters["TRANSFER"], G_NPRNG())
        l = af.obtain_value(self.parameters["LOSS"], G_NPRNG())
        i = af.obtain_value(self.parameters["INVERSION"], G_NPRNG())
        p = af.obtain_value(self.parameters["TRANSPOSITION"], G_NPRNG())
        o = af.obtain_value(self.parameters["ORIGINATION"], G_NPRNG())

        return d,t,l,i,p,o

    def generate_empirical_rates(self):

        mlen = len(self.empirical_rates)
        d,t,l = self.empirical_rates[G_NPRNG().integers(mlen)]

        return d,t,l

    def read_rates(self, rates_folder: Path):

        self.branch_event_rates = dict()
        self.branch_extension_rates = dict()
        self.transfer_rates = dict()

        with open(rates_folder / EVENTRATES) as f:
            f.readline()
            for line in f:
                sp, d, t, l, i, c, o,  = line.split("\t")
                self.branch_event_rates[sp] = tuple([float(x) for x in (d, t, l, i, c, o)])

        with open(rates_folder / EXTENSIONRATES) as f:
            f.readline()
            for line in f:
                sp, d, t, l, i, c,  = line.split("\t")
                self.branch_extension_rates[sp] = tuple([x for x in (d, t, l, i, c)])

        with open(rates_folder / TRANSFERRATES) as f:
            f.readline()
            for line in f:
                dn, rc, wt  = line.split("\t")
                if dn not in self.transfer_rates:
                    self.transfer_rates[dn] = dict()
                if rc not in self.transfer_rates[dn]:
                    self.transfer_rates[dn][rc] = 0.0

                self.transfer_rates[dn][rc] = float(wt)


    def choose_event(self, tandemdup, duplication, transfer, loss, inversion,
                     transposition, origination):

        draw = G_NPRNG().choice([TDUP, DUP, FER, LOSS, INV, POS, ORIG], 1,
                                p=af.normalize([tandemdup, duplication, transfer, loss,
                                                inversion, transposition, origination]))
        return draw


    def choose_event_i(self, duplication, transfer, loss, inversion,
                       transposition, origination, remove, rewire):

        draw = G_NPRNG().choice([TDUP, FER, LOSS, INV, POS, ORIG, "RM", "RW"], 1,
                                p=af.normalize([duplication, transfer, loss,
                                                inversion, transposition,
                                                origination, remove, rewire]))
        return draw


    def choose_recipient(self, lineages_alive, donor):
        possible_recipients = sorted(x for x in lineages_alive if x != donor)
        if len(possible_recipients) > 1:
            recipient = G_RNG().choice(possible_recipients)
            return recipient
        else:
            return None


    def evolve_genomes(self, tandemdup, duplication, transfer, loss, inversion,
                       transposition, origination,
                       time: float) -> None | tuple[str, str | tuple[str, str]]:
        """
        Choose a random rearrangement and apply it to one of the active
        lineages.

        Returns
        -------
        None | Tuple[str, str]
            (type, lineage), where `lineage` is the species of the genome
            where it was applied.  None if no rearrangement could be applied.
        """
        lineage = G_RNG().choice(sorted(self.active_genomes))

        event = self.choose_event(tandemdup, duplication, transfer, loss,
                                  inversion, transposition, origination)

        if event == TDUP:
            d_e = self.parameters["TANDEMDUP_EXTENSION"]
            self.make_tandemdup(d_e, lineage, time)
            return TDUP, lineage

        elif event == DUP:
            u_e = self.parameters["DUPLICATION_EXTENSION"]
            self.make_duplication(u_e, lineage, time)
            return DUP, lineage

        elif event == FER:
            t_e = self.parameters["TRANSFER_EXTENSION"]

            possible_recipients = sorted(x for x in self.active_genomes
                                         if x != lineage)

            if len(possible_recipients) > 0:

                donor = lineage

                # We choose a recipient

                if self.parameters["ASSORTATIVE_TRANSFER"]:
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = G_RNG().choice(possible_recipients)

                self.make_transfer(t_e, donor, recipient, time)
                return FER, (donor, recipient)

            else:
                return None

        elif event == LOSS:
            l_e = self.parameters["LOSS_EXTENSION"]

            self.make_loss(l_e, lineage, time)
            return LOSS, lineage

        elif event == INV:
            i_e = self.parameters["INVERSION_EXTENSION"]
            self.make_inversion(i_e, lineage, time)
            return INV, lineage

        elif event == POS:
            c_e = self.parameters["TRANSPOSITION_EXTENSION"]
            self.make_transposition(c_e, lineage, time)
            return POS, lineage

        elif event == ORIG:
            gene, _ = self.make_origination(lineage, time)
            chromosome = self.node_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()

            event = EventOneBreakpoint(position, 1, chromosome, ORIG, lineage, time)
            #chromosome.geneorder_history.append(event)
            self.geneorder_events[lineage].append(event)

            chromosome.insert_segment(position, [gene])
            return ORIG, lineage


    def evolve_genomes_i(self, duplication, transfer, loss, inversion,
                         transposition, origination, remove, rewire,
                         time: float) -> None | tuple[str, str | tuple[str, str]]:
        """
        Evolve genome with interactome-level rearrangements.
        """
        d_e = self.parameters["DUPLICATION_EXTENSION"]
        t_e = self.parameters["TRANSFER_EXTENSION"]
        l_e = self.parameters["LOSS_EXTENSION"]
        i_e = self.parameters["INVERSION_EXTENSION"]
        c_e = self.parameters["TRANSPOSITION_EXTENSION"]

        lineage = G_RNG().choice(sorted(self.active_genomes))

        event = self.choose_event_i(duplication, transfer, loss, inversion,
                                    transposition, origination, remove, rewire)

        if event == TDUP:

            self.make_tandemdup_interactome(d_e, lineage, time)
            return TDUP, lineage

        elif event == FER:

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                donor = lineage
                if self.parameters["ASSORTATIVE_TRANSFER"]:
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = G_RNG().choice(possible_recipients)

                self.make_transfer_interactome(t_e, donor, recipient, time)
                return FER, (donor, recipient)

            else:
                return None

        elif event == LOSS:

            self.make_loss_interactome(l_e, lineage, time)
            return LOSS, lineage

        elif event == INV:
            self.make_inversion(i_e, lineage, time)
            return INV, lineage

        elif event == POS:
            self.make_transposition(c_e, lineage, time)
            return POS, lineage

        elif event == ORIG:

            gene, _ = self.make_origination(lineage, time)
            chromosome = self.node_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            chromosome.insert_segment(position, [gene])

            event = EventOneBreakpoint(position, 1, chromosome, ORIG, lineage, time)
            #chromosome.geneorder_history.append(event)
            self.geneorder_events[lineage].append(event)

            # We need to insert the gene in the interactome too, with preferential attachment

            interactome = self.node_genomes[lineage].interactome
            assert interactome

            node_degrees = [d + 1 for n, d in interactome.degree()]    #type: ignore
            choice = G_NPRNG().choice(sorted(interactome.nodes), 1,      #type: ignore
                                      p=af.normalize(node_degrees))[0]

            interactome.add_node(str(gene))
            interactome.add_edge(str(gene), choice)

            return ORIG, lineage

        elif event == "RM":

            self.make_remove_edge(lineage, time)
            return "RM", lineage


        elif event == "RW":

            self.make_rewiring_edge(lineage, time)

            return "RW", lineage


    def evolve_genomes_m(self, time) -> None | tuple[str, str | tuple[str, str]]:
        """
        Evolve genome with family-wise rates.
        """
        d_e = self.parameters["TANDEMDUP_EXTENSION"]
        u_e = self.parameters["DUPLICATION_EXTENSION"]
        t_e = self.parameters["TRANSFER_EXTENSION"]
        l_e = self.parameters["LOSS_EXTENSION"]
        i_e = self.parameters["INVERSION_EXTENSION"]
        c_e = self.parameters["TRANSPOSITION_EXTENSION"]

        ####


        ####

        mactive_genomes = sorted(self.active_genomes)
        mweights = list()
        for genome in mactive_genomes:
            lineage_weight = 0
            for chromosome in self.node_genomes[genome]:
                for gene in chromosome:
                    for r,vl in self.all_gene_families[gene.family].rates.items():
                        lineage_weight += vl
            mweights.append(lineage_weight)

        lineage = G_NPRNG().choice(mactive_genomes, 1, p=af.normalize(mweights))[0]

        d, u, t, l, i, p, o = 0, 0, 0, 0, 0, 0, 0

        for chromosome in self.node_genomes[lineage]:
            for gene in chromosome:
                d += self.all_gene_families[gene.family].rates["TANDEMDUP"]
                #u += self.all_gene_families[gene.family].rates["DUPLICATION"]
                t += self.all_gene_families[gene.family].rates["TRANSFER"]
                l += self.all_gene_families[gene.family].rates["LOSS"]

        i += af.obtain_value((self.parameters["INVERSION"]), G_NPRNG())
        p += af.obtain_value((self.parameters["TRANSPOSITION"]), G_NPRNG())
        o += af.obtain_value((self.parameters["ORIGINATION"]), G_NPRNG())

        #print(d,t,l,i,p,o)

        event = self.choose_event(d, u, t, l, i, p, o)

        ####

        if event == TDUP:

            self.make_tandemdup(d_e, lineage, time, family_mode=True)
            return TDUP, lineage

        elif event == DUP:

            self.make_duplication(u_e, lineage, time, family_mode=True)
            return DUP, lineage


        elif event == FER:

            # We choose a recipient

            possible_recipients = sorted(x for x in self.active_genomes
                                         if x != lineage)

            if len(possible_recipients) > 0:

                donor = lineage
                if self.parameters["ASSORTATIVE_TRANSFER"]:
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = G_RNG().choice(possible_recipients)

                self.make_transfer(t_e, donor, recipient, time, family_mode=True)
                return FER, (donor, recipient)

            else:
                return None

        elif event == LOSS:
            self.make_loss(l_e, lineage, time, family_mode=True)
            return LOSS, lineage

        elif event == INV:
            self.make_inversion(i_e, lineage, time)
            return INV, lineage

        elif event == POS:
            self.make_transposition(c_e, lineage, time)
            return POS, lineage

        elif event == ORIG:

            if  self.parameters["RATE_FILE"] == "False":
                gene, _ = self.make_origination(lineage, time, family_mode=True)

            else:
                gene, _ = self.make_origination(lineage, time, family_mode=True,
                                                empirical_rates=True)

            chromosome = self.node_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            chromosome.insert_segment(position, [gene])

            event = EventOneBreakpoint(position, 1, chromosome, ORIG, lineage, time)
            #chromosome.geneorder_history.append(event)
            self.geneorder_events[lineage].append(event)

            return ORIG, lineage


    def advanced_evolve_genomes(self, time) -> None | tuple[str, str | tuple[str, str]]:
        """
        Evolve genome with user-defined rates (using RateCustomizer).
        """
        active_genomes = sorted(self.active_genomes)
        lineage = G_NPRNG().choice(active_genomes, 1, p=af.normalize(
            [sum(self.branch_event_rates[x]) for x in active_genomes]))[0]

        d,t,l,i,c,o = self.branch_event_rates[lineage]
        u = af.obtain_value('f:0', G_NPRNG())  #TODO: add functionality to RateCustomizer

        event = self.choose_event(d,u,t,l,i,c,o)

        d_e, t_e, l_e, i_e, c_e = self.branch_extension_rates[lineage]
        u_e = 'f:0'         #TODO: add functionality to RateCustomizer

        if event == TDUP:
            self.make_tandemdup(d_e, lineage, time)
            return TDUP, lineage

        elif event == FER:

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                recipient = self.choose_advanced_recipient(possible_recipients, lineage)
                if recipient != None:
                    donor = lineage
                    self.make_transfer(t_e, donor, recipient, time)
                    return FER, (donor, recipient)

            else:
                return None

        elif event == LOSS:

            self.make_loss(l_e, lineage, time)
            return LOSS, lineage

        elif event == INV:
            self.make_inversion(i_e, lineage, time)
            return INV, lineage

        elif event == POS:
            self.make_transposition(c_e, lineage, time)
            return POS, lineage

        elif event == ORIG:

            gene, _ = self.make_origination(lineage, time)

            chromosome = self.node_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            chromosome.insert_segment(position, [gene])

            event = EventOneBreakpoint(position, 1, chromosome, ORIG, lineage, time)
            #chromosome.geneorder_history.append(event)
            self.geneorder_events[lineage].append(event)

            return ORIG, lineage


    def advanced_evolve_genomes_f(self, tandemdup, duplication, transfer, loss,
                                  inversion, transposition, origination,
                                  time) -> None | tuple[str, str | tuple[str, str]]:
        """
        Evolve genome with intergenes.
        """
        # Evolve genomes with intergenes
        d_e = int(af.obtain_value(self.parameters["TANDEMDUP_EXTENSION"], G_NPRNG()))
        u_e = int(af.obtain_value(self.parameters["DUPLICATION_EXTENSION"], G_NPRNG()))
        t_e = int(af.obtain_value(self.parameters["TRANSFER_EXTENSION"], G_NPRNG()))
        l_e = int(af.obtain_value(self.parameters["LOSS_EXTENSION"], G_NPRNG()))
        i_e = int(af.obtain_value(self.parameters["INVERSION_EXTENSION"], G_NPRNG()))
        c_e = int(af.obtain_value(self.parameters["TRANSPOSITION_EXTENSION"], G_NPRNG()))

        distribution  = self.parameters["GENE_LENGTH"].split(":")[0]

        if distribution == "f" or distribution == "n":
            mean_gene_length = int(self.parameters["GENE_LENGTH"].split(":")[1].split(";")[0])
        elif distribution == "u":
            u1,u0 = self.parameters["GENE_LENGTH"].split(":")[1].split(";")
            mean_gene_length = (int(u1) - int(u0))/2
        else:
            print("Error, please switch the distribution type for the gene length")
            return None

        mean_intergene_length = int(self.parameters["INTERGENE_LENGTH"])
        multiplier = 1.0 / mean_intergene_length

        lineage = G_RNG().choice(sorted(self.active_genomes))

        event = self.choose_event(tandemdup, duplication, transfer, loss,
                                  inversion, transposition, origination)

        self.update_genome_indices(lineage)

        try: 
            r = self.select_advanced_length(lineage, 1/d_e * multiplier)
        except CoordinateChoiceError:
            return None
                

        ch, c1, c2, d = r
        d = T_DIR.RIGHT

        if event == TDUP:

            if self.make_tandemdup_intergenic(ch, c1, c2, d, lineage, time):
                return TDUP, lineage
            else:
                return None

        elif event == FER:

            # We choose a recipient

            possible_recipients = sorted(x for x in self.active_genomes
                                         if x != lineage)

            if len(possible_recipients) > 0:
                donor = lineage
                if self.parameters["ASSORTATIVE_TRANSFER"]:
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = G_RNG().choice(possible_recipients)

                try: 
                    r = self.select_advanced_length(lineage, 1/t_e * multiplier)
                except CoordinateChoiceError:
                    return None

                ch, c1, c2, d = r

                chreceptor = self.node_genomes[recipient].select_random_chromosome()
                assert isinstance(chreceptor, CircularChromosome)
                chreceptor.update_locations()
                c3 = chreceptor.select_random_coordinate_in_intergenic_regions()

                if self.make_transfer_intergenic(ch, c1, c2, d, donor, chreceptor,
                                                 c3, recipient, time):
                    return FER, (donor, recipient)
                else:
                    return None

            else:
                return None

        elif event == LOSS:

            pseudo = False
            if G_NPRNG().uniform(0,1) <= float(self.parameters["PSEUDOGENIZATION"]):
                pseudo = True

            if self.make_loss_intergenic(ch, c1, c2, d, lineage, time, pseudo):
                return LOSS, lineage
            else:
                return None

        elif event == INV:

            if self.make_inversion_intergenic(ch, c1, c2, d, lineage, time):
                return INV, lineage
            else:
                return None

        elif event == POS:

            ch, c1, c2, d = r
            c3 = ch.select_random_intergenic_coordinate_excluding(c1, c2, d)

            if c3 == None:
                return None

            if self.make_transposition_intergenic(ch, c1, c2, d, c3, lineage, time):
                return POS, lineage
            else:
                return None

        elif event == ORIG:
                
            ch = self.node_genomes[lineage].select_random_chromosome()
            intergene_coordinate = ch.select_random_coordinate_in_intergenic_regions()
            self.make_origination_intergenic(ch, intergene_coordinate,
                                             lineage, time)
            return ORIG, lineage


#    def advanced_evolve_genomes_f_debug(self, duplication, transfer, loss, inversion, transposition, origination, time):
#        
#        
#
#        for ch in self.node_genomes[lineage]:
#            ch.obtain_flankings()
#            ch.obtain_locations()
#
#        if event == TDUP:
#
#            r = self.select_advanced_length(lineage, 1/d_e * multiplier)
#            if r == None:
#                return None
#            else:
#                
#                #ch, c1, c2, d = r
#                c1, c2, d = r
#                self.make_duplication_within_intergene(ch, c1, c2, d, lineage, time)
#
#            return TDUP, lineage
#
#        elif event == FER:
#
#
#            # We choose a recipient
#
#            possible_recipients = sorted(x for x in self.active_genomes if x != lineage)
#
#            if len(possible_recipients) > 0:
#
#                if self.parameters["ASSORTATIVE_TRANSFER"]:
#                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
#                    if recipient == None:
#                        return None
#                else:
#                    recipient = random.choice(possible_recipients)
#
#                donor = lineage
#
#                r = self.select_advanced_length(lineage, 1/t_e * multiplier)
#
#                if r == None:
#                    return None
#                else:
#                    c1, c2, d = r
#                    self.make_transfer_intergenic(ch, c1, c2, d, donor, recipient, time)
#
#                return FER, donor + "->" + recipient
#
#            else:
#                return None
#
#        elif event == LOSS:
#
#            r = self.select_advanced_length(lineage, 1/l_e * multiplier)
#
#            if r == None:
#                return None
#            else:
#                c1, c2, d = r
#                pseudo = False
#                if numpy.random.uniform(0,1) <= float(self.parameters["PSEUDOGENIZATION"]):
#                    pseudo = True
#                self.make_loss_intergenic(ch, c1, c2, d, lineage, time, pseudo)
#
#            return LOSS, lineage
#
#        elif event == INV:
#
#            r = self.select_advanced_length(lineage, 1/i_e * multiplier)
#
#            if r == None:
#                return None
#            else:
#                c1, c2, d = r
#                self.make_inversion_intergenic(ch, c1, c2, d, lineage, time)
#
#            return INV, lineage
#
#        elif event == POS:
#
#            r = self.select_advanced_length(lineage, 1/c_e * multiplier)
#            if r == None:
#                return None
#
#            c1, c2, d = r
#
#            c3 = ch.select_random_intergenic_coordinate_excluding(c1, c2, d)
#            self.make_transposition_intergenic(ch, c1, c2, d, c3, lineage, time)
#
#            return POS, lineage
#
#        elif event == ORIG:
#                
#            ch = self.node_genomes[lineage].select_random_chromosome()
#            intergene_coordinate = ch.select_random_coordinate_in_intergenic_regions()
#            self.make_origination_intergenic(ch, intergene_coordinate,
#                                             lineage, time)
#
#            return ORIG, lineage


    def get_time_to_next_event(self, n: int, events: list[float]):

        total = 0.0
        for __ in range(n):
            total += sum(events)

        if total == 0:
            return 1000000000000000.0 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        else:
            time = G_NPRNG().exponential(1 / total)
            return time

    def get_time_to_next_event_advanced_modes(self):
        # To obtain the time to next event in case that we have different rates per branch
        total = 0.0

        for lineage in self.active_genomes:
            total += sum(self.branch_event_rates[lineage])

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        time = G_NPRNG().exponential(1 / total)
        return time

    def get_time_to_next_event_family_mode(self):

        total = 0.0

        for lineage in self.active_genomes:
            for chromosome in self.node_genomes[lineage]:
                for gene in chromosome:
                    for r,vl in self.all_gene_families[gene.family].rates.items():
                        total += vl

        total_active =  len(self.active_genomes)

        total += af.obtain_value((self.parameters["INVERSION"]), G_NPRNG()) * total_active
        total += af.obtain_value((self.parameters["TRANSPOSITION"]), G_NPRNG()) * total_active
        total += af.obtain_value((self.parameters["ORIGINATION"]), G_NPRNG()) * total_active

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        time = G_NPRNG().exponential(1 / total)

        return time


    def increase_distances(self, time_to_next_event, active_lineages):

        for node in active_lineages:
            node.dist += time_to_next_event


    def make_origination(self, species_tree_node, time, family_mode = False,
                         empirical_rates = False) -> tuple[Gene, GeneFamily]:

        self.gene_families_counter += 1
        gene_family_id = str(self.gene_families_counter)

        gene = Gene()
        gene.determine_orientation()

        gene.family = str(self.gene_families_counter)
        gene.species = species_tree_node

        gene_family = GeneFamily(gene_family_id, time)
        gene_family.length = int(af.obtain_value(self.parameters["GENE_LENGTH"], G_NPRNG()))
        gene_family.initial_orientation = gene.orientation
        
        gene.length = gene_family.length

        gene_family.genes.append(gene)
        gene.gene_id = gene_family.obtain_new_gene_id()

        self.all_gene_families[gene_family_id] = gene_family
        self.all_gene_families[gene.family].append_event(time, ORIG, species_tree_node)

        if family_mode == True and empirical_rates == False:

            d, t,l, _, _, _ = self.generate_new_rates()
            gene_family.rates["TANDEMDUP"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        elif family_mode == True and empirical_rates == True:
            d, t, l = self.generate_empirical_rates()
            gene_family.rates["TANDEMDUP"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        return gene, gene_family


    def make_origination_intergenic(self, chromosome: Chromosome, c: int,
                                    lineage, time) -> Gene:

        if isinstance(chromosome, LinearChromosome):
            raise NotImplementedError('Origination in intergenes not implemented'
                                      ' for linear chromosomes yet.')
        assert isinstance(chromosome, CircularChromosome)
        
        gene, _ = self.make_origination(lineage, time)
        
        location, gene_position = chromosome.get_location_from_coord(c, True)
        chromosome.insert_gene_within_intergene(c, location, gene)
 
        #Nucleotide-level origination event
        orig = Origination(location, c, gene.length, gene.family, gene.orientation, lineage, time)
        chromosome.event_history.append(orig)

        #Gene-order level origination event
        event = EventOneBreakpoint(gene_position, 1, chromosome, ORIG, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)

        return gene


    def make_speciation(self, sp:str, c1:str, c2:str, time:float,
                        intergene=False):
        """
        Make a speciation event at species node `sp`, generating two new
        genomes for the child species `c1` and `c2`.
        """
        if self.save_all:
            self.add_genome_to_all_genomes(sp, True)

        if sp not in self.geneorder_events:
            self.geneorder_events[sp] = []

        genome_sp = self.node_genomes[sp]

        genome1 = Genome()
        genome2 = Genome()

        if hasattr(genome_sp, 'interactome'):
            genome1.interactome = copy.deepcopy(genome_sp.interactome)
            genome2.interactome = copy.deepcopy(genome_sp.interactome)

            new_names_1 = dict()
            new_names_2 = dict()

        for chromosome in genome_sp:

            shape = chromosome.shape

            if shape == "C":
                ch1 = CircularChromosome(chromosome.name)
                ch2 = CircularChromosome(chromosome.name)

            else:
                assert shape == "L"
                raise(NotImplementedError("Linear chromosomes not implemented"))
                ch1 = LinearChromosome(chromosome.name)
                ch2 = LinearChromosome(chromosome.name)

            genome1.chromosomes.append(ch1)
            genome2.chromosomes.append(ch2)

            if chromosome.has_intergenes:
                ch1.has_intergenes = True
                ch2.has_intergenes = True

            for gene in chromosome:

                new_gene1 = self.copy_segment([gene])[0]
                new_gene2 = self.copy_segment([gene])[0]

                new_gene1.species = c1
                new_gene2.species = c2

                new_gene1.orientation = gene.orientation
                new_gene2.orientation = gene.orientation

                new_gene1.length = gene.length
                new_gene2.length = gene.length

                gene_family = self.all_gene_families[gene.family]
                gene_family.genes.append(new_gene1)
                gene_family.genes.append(new_gene2)

                new_gene1.family = gene.family
                new_gene2.family = gene.family

                ch1.genes.append(new_gene1)
                ch2.genes.append(new_gene2)

                if hasattr(genome_sp, 'interactome'):

                    new_names_1[str(gene)] = str(new_gene1)
                    new_names_2[str(gene)] = str(new_gene2)

                gene.active = False

                # The code for the node is:
                # 1. Branch of the species tree that splits
                # 2. Id of the gene that is split
                # 3. Branch of the species tree first child
                # 4. Id of the gene that goes to first child
                # 5. Branch of the species tree second child
                # 6. Id of the gene that goes to second child

                nodes = [sp,
                         gene.gene_id,
                         c1,
                         new_gene1.gene_id,
                         c2,
                         new_gene2.gene_id
                         ]

                self.all_gene_families[gene.family].append_event(time, "S", ";".join(map(str,nodes)))

            for intergene in chromosome.intergenes:

                new_intergene1 = Intergene()
                new_intergene2 = Intergene()
                new_intergene1.length = intergene.length
                new_intergene2.length = intergene.length
                ch1.intergenes.append(new_intergene1)
                ch2.intergenes.append(new_intergene2)

        genome1.update_genome_species(c1)
        genome2.update_genome_species(c2)

        if hasattr(genome_sp, 'interactome'):
            genome1.interactome = nx.relabel_nodes(genome1.interactome, new_names_1)
            genome2.interactome = nx.relabel_nodes(genome2.interactome, new_names_2)
            #nx.relabel_nodes(genome1.interactome, new_names_1)
            #nx.relabel_nodes(genome2.interactome, new_names_2)

        return genome1, genome2


    def make_extinction(self, sp:str, time:float):
        """
        Deactivate the genome at species node `sp` due to extinction.
        """
        if self.save_all:
            self.add_genome_to_all_genomes(sp, True)

        if sp not in self.geneorder_events:
            self.geneorder_events[sp] = []

        genome = self.node_genomes[sp]

        for chromosome in genome:
            for gene in chromosome:
                gene.active = False
                self.all_gene_families[gene.family].append_event(time, "E", ";".join(map(str,[sp, gene.gene_id])))


    def make_end(self, time):
        """
        Clean things up for any remaining active genomes at the end of the simulation.
        """
        for genome_name in self.active_genomes:
            genome = self.node_genomes[genome_name]
            for chromosome in genome:
                for gene in chromosome:
                    gene.active = False
                    self.all_gene_families[gene.family].append_event(time, "F", ";".join(
                        map(str, [genome.species, gene.gene_id])))

            if self.save_all:
                self.add_genome_to_all_genomes(genome_name, True)

            if genome_name not in self.geneorder_events:
                self.geneorder_events[genome_name] = []


    def make_tandemdup(self, p: float, lineage: str, time: float,
                       family_mode=False):
        """
        Tandemly duplicate a segment of genes.  This is always done to the right
        of the original segment.
        """
        self.make_tandemdup_core(p, lineage, time, family_mode)


    def make_tandemdup_core(self, p: float, lineage, time, family_mode=False) \
        -> tuple[list[Gene], list[Gene], list[Gene]]:
        """
        Tandemly duplicate a segment of genes on the current lineage (i.e. on
        'self.node_genomes').  This is always done to the right of the original
        segment.

        Returns
        -------
        tuple[list[Gene], list[Gene], list[Gene]]
            (original segment, left_copy, right_copy)
            Where original segment is the part that was duplicated, and replaced
            with left_copy + right_copy.
        """
        chromosome = self.node_genomes[lineage].select_random_chromosome()

        if family_mode:
            affected_indices = chromosome.obtain_affected_indices_family_rates(p, self.all_gene_families, "TANDEMDUP")
        else:
            affected_indices = chromosome.obtain_affected_indices(p)

        event = EventTwoBreakpoints(affected_indices, chromosome, TDUP, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)
        segment = chromosome.obtain_segment(affected_indices)

        # Now we create two segments

        copied_segment1 = self.copy_segment(segment)
        copied_segment2 = self.copy_segment(segment)

        # Replace the original segment with the first of the copies
        chromosome.replace_segment(affected_indices, copied_segment1)

        # Insert the two new segments after the last position of the old segment
        chromosome.insert_segment(affected_indices[-1]+1, copied_segment2)

        # We have to register in the affected gene families that there has been a duplication

        for i, gene in enumerate(segment):

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            gene.active = False

            # We add the genes to the list of genes in the gene family
            gene_family = gene.family

            self.all_gene_families[gene_family].genes.append(copied_segment1[i])
            self.all_gene_families[gene_family].genes.append(copied_segment2[i])
            self.all_gene_families[gene_family].append_event(time, TDUP, ";".join(map(str, nodes)))

        return segment, copied_segment1, copied_segment2


    def make_duplication(self, p: float, lineage: str, time: float,
                         family_mode = False):
        """
        Duplicate a segment of genes to a (uniform) random location.
        """
        chromosome = self.node_genomes[lineage].select_random_chromosome()
        assert isinstance(chromosome, CircularChromosome)

        if family_mode == True:
            raise(NotImplementedError("Family mode not implemented for duplication"))
            affected_indices = chromosome.obtain_affected_genes_accounting_for_family_rates(p, self.all_gene_families, "TANDEMDUP")
        else:
            affected_indices = chromosome.obtain_affected_indices(p)
            while len(affected_indices) == len(chromosome.genes):
                affected_indices = chromosome.obtain_affected_indices(p)

        segment = chromosome.obtain_segment(affected_indices)

        # Now we create two segments
        copied_segment1 = self.copy_segment(segment)
        copied_segment2 = self.copy_segment(segment)

        # Replace the segment with the first copy
        chromosome.replace_segment(affected_indices, copied_segment1)

        #if([g.family for g in segment] != \
        #   [chromosome.genes[i].family for i in affected_indices]):
        #    print(affected_indices[0], affected_indices[-1])
        #    print(f'{segment} != {chromosome.genes[affected_indices[0]:affected_indices[-1]+1]}')

        # Insert the second segment at a random location.
        insert_pos = chromosome.select_random_position(affected_indices)
        chromosome.insert_segment(insert_pos, copied_segment2)

        event = EventThreeBreakpoints(affected_indices, insert_pos, chromosome,
                                      DUP, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)

        # Register the duplication in the affected gene families
        for i, gene in enumerate(segment):
            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            gene.active = False

            # We add the genes to the list of genes in the gene family
            gene_family = gene.family

            self.all_gene_families[gene_family].genes.append(copied_segment1[i])
            self.all_gene_families[gene_family].genes.append(copied_segment2[i])
            self.all_gene_families[gene_family].append_event(time, DUP, ";".join(map(str, nodes)))


    def make_tandemdup_interactome(self, p, lineage, time, family_mode=False):

        segment, left_segment, right_segment = self.make_tandemdup_core(p, lineage, time, family_mode)

        # We need to create two dicts to change the names of the interactome
        new_genes_1 = dict()

        interactome = self.node_genomes[lineage].interactome
        for i, gene in enumerate(segment):
            # We update the new gene name (which by default is going to be the copied_segment1)

            gene1 = left_segment[i]
            gene2 = right_segment[i]

            new_genes_1[str(gene)] = str(gene1)

            interactome = nx.relabel_nodes(interactome, new_genes_1)

            # WE ADD THE NEW NODE

            interactome.add_node(str(gene2))

            # We distribute node depending on the parameter PROPORTION

            # If p is 1, all the links go to the first node
            # If p is 0, all the links go to the second node
            # If p is 0.5, equal repartition.

            PROPORTION = 0.5

            n_edges_to_old_node = int(PROPORTION * len(interactome.edges(str(gene1))))
            myedges = list(interactome.edges(str(gene1)))
            G_RNG().shuffle(myedges)
            edges_to_new_node = myedges[n_edges_to_old_node:]
            edges_to_add_to_new_node = [(str(gene2), x[1]) for x in edges_to_new_node]
            edges_to_remove_to_old_node = edges_to_new_node

            # Now, I have to remove the edges

            interactome.remove_edges_from(edges_to_remove_to_old_node)
            interactome.add_edges_from(edges_to_add_to_new_node)


    def make_tandemdup_intergenic(self, chromosome: Chromosome,
                                  c1: int, c2: int, d: T_DIR,
                                  lineage: str, time: float):
        """
        Do a duplication that acts on the given pair of intergene specific
        breakpoint coordinates. Consider intergene I and J such that c1 lands in
        intergene I and c2 lands in intergene J, with gene/intergene segment S
        between the two. Then we have sequence

            I S J

        where I is split at `c1` into I0 I1 and J is split at `c2` into J0 J1.
        Then we get

            I0 I1 S J0 J1

        and the tandem duplication produces

            I0 I1 S J0 I1 S J0 J1.

        Parameters
        ----------
        c1 : int
            the first intergene specific breakpoint coordinate
        c2 : int
            the second intergene specific breakpoint coordinate
        d : T_DIR
            the direction, either RIGHT or LEFT
        lineage : str
            the lineage, which is the name of the pendant node
        time : float
            the time stamp of the event
        """
        if isinstance(chromosome, LinearChromosome):
            raise NotImplementedError('Tandem duplication within intergenes not'
                                      ' implemented for linear chromosomes.')

        try:
            r = chromosome.get_affected_region(c1, c2, d)
        except CoordinateChoiceError:
            return False


        genepositions, intergenepositions, leftlengths, rightlengths, int1, int2 = r
        segment = chromosome.obtain_segment(genepositions)
        #intergene_segment = chromosome.obtain_intergenic_segment(intergenepositions[1:])

        # We duplicate the genes

        new_segment_1 = self.copy_segment(segment)
        new_segment_2 = self.copy_segment(segment)

        # And the intergenes, which will have incorrect lengths until we update
        # them later

        new_intergene_segment_1 = [copy.deepcopy(chromosome.intergenes[x])
                                   for x in intergenepositions[1:]]
        new_intergene_segment_2 = [copy.deepcopy(chromosome.intergenes[x])
                                   for x in intergenepositions[1:]]

        scar0 = chromosome.intergenes[intergenepositions[0]]
        scar1 = new_intergene_segment_1[-1]
        scar2 = new_intergene_segment_2[-1]

        ###
        ###

            # Get old lengths from last intergene before modifying chromosome.
        specificlen = chromosome.intergenes[-1].specific_flanking[1]
        totallen = chromosome.intergenes[-1].total_flanking[1]


            #Replace the original set of genes/intergenes with the first copies:
        for pos, newgene in zip(genepositions, new_segment_1):
            chromosome.genes[pos] = newgene
        for pos, newintergene in zip(intergenepositions[1:], new_intergene_segment_1):
            chromosome.intergenes[pos] = newintergene
            
        position = genepositions[-1] + 1

            #Insert the second copy of new genes and intergenes after the first.
            #map_of_locations will later be updated by obtain_locations().
        for i, gene in enumerate(new_segment_2):
            chromosome.genes.insert(position + i, gene)
        for i, intergene in enumerate(new_intergene_segment_2):
            chromosome.intergenes.insert(position + i, intergene)

        # We adjust the new intergenes lengths and record the TandemDup:

        if d == T_DIR.LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        #Nucleotide-level tandem duplication event
        dup = TandemDup(int1, int2, c1, c2, len(new_intergene_segment_1),
                        specificlen, totallen, lineage, time)

        #Gene-order level tandem duplication event
        event = EventTwoBreakpoints(genepositions, chromosome, TDUP, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)
            
        scar1.length = leftlengths[1] + rightlengths[0]
        scar2.length = rightlengths[1] + rightlengths[0]

        assert scar1.length == len(dup.afterC)
        assert scar2.length == len(dup.afterR)
        chromosome.event_history.append(dup)
        self._dupAssert(dup, scar0, scar1, scar2, chromosome)   #TODO: Temporary

        for i, gene in enumerate(segment):
            nodes = [gene.species,
                     gene.gene_id,
                     new_segment_1[i].species,
                     new_segment_1[i].gene_id,
                     new_segment_2[i].species,
                     new_segment_2[i].gene_id]

            gene.active = False

            gene_family = gene.family

            self.all_gene_families[gene_family].genes.append(new_segment_1[i])
            self.all_gene_families[gene_family].genes.append(new_segment_2[i])
            self.all_gene_families[gene_family].append_event(time, TDUP, ";".join(map(str, nodes)))

        return True


    def _dupAssert(self, dup: TandemDup, ileft: Intergene, icenter: Intergene,
                   iright: Intergene, chromosome: Chromosome):
        """
        Do sanity checks on the TandemDup by comparing it to the given center
        and right intergenes on the give `chromosome`.

        Parameters
        ----------
        dup : TandemDup
            the duplication
        ileft : Intergene
            the new left breakpoint from the chromosome
        icenter : Intergene
            the new center breakpoint from the chromosome
        iright : Intergene
            the new right breakpoint from the chromosome
        chromosome : Chromosome
            the chromosome that was modified
        """
        first, second = ileft.specific_flanking
        chromosome.update_flankings()

            #Specific coordinate asserts:
        assert ileft.specific_flanking[0] == dup.afterL.sc1, \
               f'{ileft.specific_flanking[0]} != {dup.afterL.sc1} ' + \
               f'{ileft.specific_flanking} {icenter.specific_flanking}'
        assert ileft.specific_flanking[1] == dup.afterL.sc2, \
               f'{ileft.specific_flanking[1]} != {dup.afterL.sc2} ' + \
               f'{icenter.specific_flanking} {iright.specific_flanking}'
        assert icenter.specific_flanking[0] == dup.afterC.sc1, \
               f'{icenter.specific_flanking[0]} != {dup.afterC.sc1} ' + \
               f'{icenter.specific_flanking} {iright.specific_flanking}'
        assert icenter.specific_flanking[1] == dup.afterC.sc2, \
               f'{icenter.specific_flanking[1]} != {dup.afterC.sc2} ' + \
               f'{icenter.specific_flanking} {iright.specific_flanking}'
        assert iright.specific_flanking[0] == dup.afterR.sc1, \
               f'{iright.specific_flanking[0]} != {dup.afterR.sc1}' + \
               f'{icenter.specific_flanking} {iright.specific_flanking}'
        assert iright.specific_flanking[1] == dup.afterR.sc2, \
               f'{iright.specific_flanking[1]} != {dup.afterR.sc2}' + \
               f'{icenter.specific_flanking} {iright.specific_flanking}'

            #Total coordinate asserts:
        assert ileft.total_flanking[0] == dup.afterL.tc1, \
               f'{ileft.total_flanking[0]} != {dup.afterL.tc1} ' + \
               f'{ileft.total_flanking} {icenter.total_flanking}'
        assert ileft.total_flanking[1] == dup.afterL.tc2, \
               f'{ileft.total_flanking[1]} != {dup.afterL.tc2} ' + \
               f'{icenter.total_flanking} {iright.total_flanking}'
        assert icenter.total_flanking[0] == dup.afterC.tc1, \
               f'{icenter.total_flanking[0]} != {dup.afterC.tc1} ' + \
               f'{icenter.total_flanking} {iright.total_flanking}'
        assert icenter.total_flanking[1] == dup.afterC.tc2, \
               f'{icenter.total_flanking[1]} != {dup.afterC.tc2} ' + \
               f'{icenter.total_flanking} {iright.total_flanking}'
        assert iright.total_flanking[0] == dup.afterR.tc1, \
               f'{iright.total_flanking[0]} != {dup.afterR.tc1}' + \
               f'{icenter.total_flanking} {iright.total_flanking}'
        assert iright.total_flanking[1] == dup.afterR.tc2, \
               f'{iright.total_flanking[1]} != {dup.afterR.tc2}' + \
               f'{icenter.total_flanking} {iright.total_flanking}'


    def choose_assortative_recipient(self, time, possible_recipients, donor):

        alpha = self.parameters["ALPHA"]
        weights = list()

        mdonor = self.complete_tree&donor

        for recipient in possible_recipients:
            mrecipient = self.complete_tree&recipient
            ca = self.complete_tree.get_common_ancestor(mrecipient, mdonor).name  #type: ignore
            x1 = self.distances_to_start[ca]
            td =  time - x1
            weights.append(td)

        beta = min(alpha * af.normalize(weights))
        val = (alpha * af.normalize(weights)) - beta
        pvector = af.normalize(numpy.exp(-val))

        draw = G_NPRNG().choice(sorted(possible_recipients), 1, p=pvector)[0]

        return draw

    def choose_advanced_recipient(self, possible_recipients, donor):

        weights = list()

        for recipient in possible_recipients:
            weights.append(self.transfer_rates[donor][recipient])

        if sum(weights) == 0:
            return None


        draw = G_NPRNG().choice(sorted(possible_recipients), 1, p=af.normalize(weights))[0]

        return draw


    def make_transfer(self, p, donor, recipient, time, family_mode = False):
        """ Do a transfer event from donor to recipient. """
        self.make_transfer_core(p, donor, recipient, time, family_mode)


    def make_transfer_core(self, p, donor, recipient, time, family_mode = False) \
        -> tuple[list[Gene], list[Gene], list[tuple[Gene, Gene]]]:
        """
        Implement the core functionality for a transfer event from donor to
        recipient (code reusability).

        Returns
        -------
        (segment, replaced_genes) : tuple[list[Gene], list[tuple[Gene, Gene]]]
            where `segment` is the list of genes that have been transferred and
            `replaced_genes` is a list of tuples (lost_gene, gained_gene)
            representing genes lost and gained in a replacement transfer; empty
            list otherwise
        """
        chromosome1 = self.node_genomes[donor].select_random_chromosome()

        if family_mode == True:
            affected_genes = chromosome1.obtain_affected_indices_family_rates(p, self.all_gene_families, "TRANSFER")
        else:
            affected_genes = chromosome1.obtain_affected_indices(p)

        segment = chromosome1.obtain_segment(affected_genes)

        inverted = False

        # Now we create two segments

        donor_segment = self.copy_segment(segment)
        recipient_segment = self.copy_segment(segment)

        # We insert the first segment (leaving transfer) in the same position
        # than the previous segment
        # We do this just to change the identifiers of the numbers

        chromosome1.insert_segment(affected_genes[0], donor_segment)

        # And we remove the old segment

        chromosome1.remove_segment(segment)

        # Now we insert the transfer segment in the recipient genome in one of
        # the homologous position.
        replaced_genes: list[tuple[Gene, Gene]] = []
        if G_NPRNG().uniform(0,1) <= self.parameters["REPLACEMENT_TRANSFER"]:

            possible_positions: list[tuple[str, Iterable[int], Chromosome]] = list()

            for chromosome in self.node_genomes[recipient]:
                for direction, positions in chromosome.get_homologous_position(segment):
                    possible_positions.append((direction, positions, chromosome))

            if len(possible_positions) != 0:
                direction, positions, chromosome2 = G_RNG().choice(sorted(possible_positions))

                positions = list(positions)
                assert isinstance(positions, Iterable), f"positions is {type(positions)}"
                assert af.well_behaved_indices(positions), f"positions are {positions}"
                loss = EventTwoBreakpoints(positions, chromosome2, LOSS, recipient, time)
                transfer = EventOneBreakpoint(positions[0], len(positions),
                                              chromosome2, AFER, recipient, time)
                #chromosome2.geneorder_history.append(loss)
                #chromosome2.geneorder_history.append(transfer)
                self.geneorder_events[recipient].append(loss)
                self.geneorder_events[recipient].append(transfer)

                if direction == BACKWARDS:
                    # Invert the segment and change orientation of the genes
                    inverted = True
                    recipient_segment = list(reversed(recipient_segment))

                    for gene in recipient_segment:
                        gene.change_orientation()
                else:
                    assert direction == FORWARDS

                # Replace the segment gene by gene:
                for i, position in enumerate(positions):
                    # First deactivate the gene
                    gene = chromosome2.genes[position]
                    gene.active = False
                    self.all_gene_families[gene.family].append_event(time, LOSS, ";".join(
                        map(str, [recipient, gene.gene_id])))

                    # And then replace
                    chromosome2.genes[position] = recipient_segment[i]
                    replaced_genes.append((gene, recipient_segment[i]))

            else:
                # Normal transfers
                chromosome2 = self.node_genomes[recipient].select_random_chromosome()
                position = chromosome2.select_random_position()
                chromosome2.insert_segment(position, recipient_segment)

                event = EventOneBreakpoint(position, len(recipient_segment),
                                           chromosome2, AFER, recipient, time)
                #chromosome2.geneorder_history.append(event)
                self.geneorder_events[recipient].append(event)
        else:
            # Normal transfer
            chromosome2 = self.node_genomes[recipient].select_random_chromosome()
            position = chromosome2.select_random_position()
            chromosome2.insert_segment(position, recipient_segment)

            event = EventOneBreakpoint(position, len(recipient_segment),
                                       chromosome2, AFER, recipient, time)
            #chromosome2.geneorder_history.append(event)
            self.geneorder_events[recipient].append(event)

        # Register the gene-order event in the donor genome
        eventcode = LFER
        if replaced_genes:
            eventcode = LFER_B if inverted else LFER_F
        event = EventTwoBreakpoints(affected_genes, chromosome1, eventcode,
                                    donor, time)

        #chromosome1.geneorder_history.append(event)
        self.geneorder_events[donor].append(event)

        # We have to register in the affected gene families that there has been a transfer event

        if inverted == True:
            # We invert again to store the event
            recipient_segment = list(reversed(recipient_segment))

        for i, gene in enumerate(segment):

            gene.active = False

            # The code for the node is:
            # 1. Branch of the species tree for the donor genome
            # 2. Id of the gene that is transferred
            # 3. Id of the gene that remains in the donor genome
            # 4. Branch of the species tree for the recipient genome
            # 5. Id of the new gene arriving

            donor_segment[i].species = donor
            recipient_segment[i].species = recipient

            nodes = [gene.species,
                     gene.gene_id,
                     donor_segment[i].species,
                     donor_segment[i].gene_id,
                     recipient_segment[i].species,
                     recipient_segment[i].gene_id]

            self.all_gene_families[gene.family].append_event(time, FER, ";".join(map(str,nodes)))

        return donor_segment, recipient_segment, replaced_genes
    

    def make_transfer_interactome(self, p, donor, recipient, time, family_mode=False):
        """
        Do a transfer event from donor to recipient, updating the
        interactome.
        """
        donor_segment, recipient_segment, replaced_genes = self.make_transfer_core(p, donor, recipient, time,
                                                                                   family_mode)

        donor_names = dict()
        recipient_names = dict()
        for i, gene in enumerate(recipient_segment):
            ## We update the interactome
            donor_names[str(gene)] = str(donor_segment[i])

            # First we update the interactome in the donor lineage
            self.node_genomes[donor].interactome = nx.relabel_nodes(self.node_genomes[donor].interactome, donor_names)

            # Second we update the interactome in the recipient lineage
            if replaced_genes:
                ## All links pass now to the new gene
                ## The old gene gets no links
                recipient_names = {str(n1):str(n2) for n1,n2 in replaced_genes}
                #new_names_2[str(gene)] = str(copied_segment2[i])
                self.node_genomes[recipient].interactome = nx.relabel_nodes(self.node_genomes[recipient].interactome, recipient_names)

            else:
                # It is not a replacement transfer. Preferential attachment
                node_degrees = [d + 1 for n, d in self.node_genomes[recipient].interactome.degree()]   #type: ignore
                choice = G_NPRNG().choice(sorted(self.node_genomes[recipient].interactome.nodes),      #type: ignore
                                                 1, p=af.normalize(node_degrees))[0]
                self.node_genomes[recipient].interactome.add_node(str(recipient_segment[i]))
                self.node_genomes[recipient].interactome.add_edge(str(recipient_segment[i]), choice)


    def make_transfer_intergenic(self, donorchrom: Chromosome,
                                 c1: int, c2: int, d: T_DIR, donor: str,
                                 receptorchrom: Chromosome, c3: int,
                                 receptor: str, time: float):
        """
        Do a transfer that acts on the given pair of intergene specific
        breakpoint coordinates.

        Parameters
        ----------
        donorchrom : Chromosome
            the donor chromosome
        c1 : int
            the first intergene specific breakpoint coordinate
        c2 : int
            the second intergene specific breakpoint coordinate
        d : T_DIR
            the direction, either RIGHT or LEFT, depending on whether the
            region goes left or right from c1
        donor : str
            the donor lineage
        receptorchrom : Chromosome
            the receptor chromosome
        c3 : int
            the intergene specific breakpoint coordinate in the receptor
        receptor : str
            the receptor lineage
        time : float
            the time stamp of the event
        """
        try:
            r = donorchrom.get_affected_region(c1, c2, d)
        except CoordinateChoiceError:
            return False

        gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r
        segment = donorchrom.obtain_segment(gpositions)

        # Get lengths from last intergene.
        specificlen = donorchrom.intergenes[-1].specific_flanking[1]
        totallen = donorchrom.intergenes[-1].total_flanking[1]
        numintergenes = len(donorchrom.intergenes)

        # Now we create two segments
        copied_segment1 = self.copy_segment(segment)
        copied_segment2 = self.copy_segment(segment)

        new_intergene_segment = [copy.deepcopy(donorchrom.intergenes[x])
                                 for x in igpositions[1:]]

        # We insert the first segment (leaving transfer) in the same position as the previous segment
        # We do this just to change the identifiers of the numbers

        # We insert in the same place

        for old_gene, new_gene in zip(segment, copied_segment1):
            insert = (donorchrom.genes).index(old_gene)
            (donorchrom.genes).pop(insert)
            (donorchrom.genes).insert(insert, new_gene)

        # We remove the old copies:

        # Normal transfer

        int3, gene_position = receptorchrom.get_location_from_coord(c3, True)
        position = int3.position + 1

        for i, gene in enumerate(copied_segment2):
            receptorchrom.genes.insert(position + i, gene)
        for i, intergene in enumerate(new_intergene_segment):
            receptorchrom.intergenes.insert(position + i, intergene)

        cut_position = (c3 - int3.sc1, int3.sc2 - c3)

        scar1 = receptorchrom.intergenes[int3.position]
        scar2 = receptorchrom.intergenes[position + i]

        if d == T_DIR.LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        scar1.length = leftlengths[1] + cut_position[0]
        scar2.length = rightlengths[0] + cut_position[1]

        tr = Transfer(int1, int2, c1, c2, specificlen, totallen, donor,
                      numintergenes, int3, c3, receptor, time)

        tr_d = Transfer(int1, int2, c1, c2, specificlen, totallen, donor,
                        numintergenes, int3, c3, receptor, time)

        #Nucleotide level origination event
        receptorchrom.event_history.append(tr)
        donorchrom.event_history.append(tr_d)

        #Gene-order level origination event
        event = EventTwoBreakpoints(gpositions, donorchrom, LFER, donor, time)
        #donorchrom.geneorder_history.append(event)
        self.geneorder_events[donor].append(event)
        event = EventOneBreakpoint(gene_position, 1, receptorchrom, AFER,
                                   receptor, time)
        #receptorchrom.geneorder_history.append(event)
        self.geneorder_events[receptor].append(event)

        tr.sister_event = tr_d
        tr_d.sister_event = tr

        tr.lineage = receptor
        tr_d.lineage = donor

        tr_d.etype = "N" 

        # We have to register in the affected gene families that there has been a transfer event


        for i, gene in enumerate(segment):

            gene.active = False

            # The code for the node is:
            # 1. Branch of the species tree for the donor genome
            # 2. Id of the gene that is transferred
            # 3. Id of the gene that remains in the donor genome
            # 4. Branch of the species tree for the recipient genome
            # 5. Id of the new gene arriving

            copied_segment1[i].species = donor
            copied_segment2[i].species = receptor

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            self.all_gene_families[gene.family].append_event(time, FER, ";".join(map(str, nodes)))

        return True


    def make_loss(self, p, lineage, time, family_mode = False):
        self.make_loss_core(p, lineage, time, family_mode)


    def make_loss_core(self, p, lineage, time, family_mode = False,
                       interactome=False) -> list[Gene]:
        chromosome = self.node_genomes[lineage].select_random_chromosome()
        if interactome:
            i = self.node_genomes[lineage].interactome
            affected_indices = chromosome.obtain_affected_genes_accounting_for_connectedness(p, i)
        if family_mode == True:
            affected_indices = chromosome.obtain_affected_indices_family_rates(p, self.all_gene_families, "LOSS")
        else:
            affected_indices = chromosome.obtain_affected_indices(p)
        segment = chromosome.obtain_segment(affected_indices)

        # Now we check we are not under the minimum size
        if len(chromosome) - len(affected_indices) <= self.parameters["MIN_GENOME_SIZE"]:
            return []

        event = EventTwoBreakpoints(affected_indices, chromosome, LOSS, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned
        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.family].append_event(time, LOSS, ";".join(map(str,[lineage, gene.gene_id])))

        return segment


    def make_loss_intergenic(self, chromosome: Chromosome, c1, c2, d: T_DIR,
                             lineage, time, pseudo=False):

        try:
            r = chromosome.get_affected_region(c1, c2, d)
        except CoordinateChoiceError:
            return False

        gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r

        segment = chromosome.obtain_segment(gpositions)
        intergene_segment = chromosome.obtain_intergenic_segment(igpositions[1:])

        # Before continuing, we need to verify that the event does not make the 
        # genome smaller than the minimum size allowed FIX --> This should be a parameter
        
        if len(chromosome.genes) <= len(segment):
            # The event does not occur
            return False
        
        # We need the adjustment factor if the event wraps,
        # and wether the genes are at the end of the chromosome
        # or at the beginning
        
        if d == T_DIR.RIGHT and c1 > c2 or d==T_DIR.LEFT and c1 < c2: # If the event wraps
            
            adjustment_factor = chromosome.genes[gpositions[-1] + 1].total_flanking[0]
            # The coordinate of the first gene not affected by the event
            
        else:
            adjustment_factor = None

        scar1 = chromosome.intergenes[igpositions[0]]

        # Get old lengths from last intergene before modifying chromosome.
        specificlen = chromosome.intergenes[-1].specific_flanking[1]
        totallen = chromosome.intergenes[-1].total_flanking[1]

        # Now we remove the genes

        for gene in segment:
            chromosome.genes.remove(gene)

        # Now we remove the intergenes

        for intergene in intergene_segment:
            chromosome.intergenes.remove(intergene)

        # We modify the length of the scar:

        if d == T_DIR.LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        pseudo_intergenes = []
        pseudo_genes = []

        if pseudo:
            pseudo_intergenes = intergene_segment
            pseudo_genes = segment

        # Nucleotide-level loss event
        loss = Loss(int1, int2, c1, c2, specificlen, totallen, lineage, time,
                    pseudo, pseudo_intergenes, pseudo_genes, adjustment_factor)
        chromosome.event_history.append(loss)
        # Gene-order level loss event
        event = EventTwoBreakpoints(gpositions, chromosome, LOSS, lineage, time)
        self.geneorder_events[lineage].append(event)

        if pseudo:

            # We need to add the length of the genes removed
            scar1.length = sum(leftlengths) + sum(rightlengths) \
                           + sum([x.length for x in segment]) \
                           + sum([x.length for x in intergene_segment[:-1]])
        else:

            scar1.length = leftlengths[0] + rightlengths[1]

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned

        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.family].append_event(time, LOSS, ";".join(map(str,[lineage, gene.gene_id])))

        return True


    def make_loss_interactome(self, p, lineage, time, family_mode=False):
        interactome = self.node_genomes[lineage].interactome
        segment = self.make_loss_core(p, lineage, time, family_mode, True)

        # Remove each from the connectome
        for gene in segment:
            interactome.remove_node(str(gene))


    def make_inversion(self, p, lineage, time):

        chromosome = self.node_genomes[lineage].select_random_chromosome()
        affected_indices = chromosome.obtain_affected_indices(p)
        segment = chromosome.obtain_segment(affected_indices)
        chromosome.invert_segment(affected_indices)

        event = EventTwoBreakpoints(affected_indices, chromosome, INV, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)

        for gene in segment:
            self.all_gene_families[gene.family].append_event(time, INV, ";".join(map(str,[lineage, gene.gene_id])))

    
    def make_inversion_intergenic(self, chromosome: Chromosome, c1: int,
                                  c2: int, d: T_DIR, lineage: str, time: float):
        
        """
        Do an inversion that acts on the given pair of intergene specific
        breakpoint coordinates. Consider intergene I and J such that c1 lands in
        intergene I and c2 lands in intergene J, with gene-intergene-gene segment G1-I1-G2
        between the two. Then we have sequence

            I G1-I1-G2 J

        where I is split at `c1` into I0 I1 and J is split at `c2` into J0 J1.
        Then we get

            I0 I1 G1-I1-G2 J0 J1

        and the inversion produces

            I0 I1 G2-I1-G1 J0 J1.

        Notes
        -----
            Must update coordinates with `self.update_genome_indices` before
            calling this.

        Parameters
        ----------
        c1 : int
            the first intergene specific breakpoint coordinate
        c2 : int
            the second intergene specific breakpoint coordinate
        d : T_DIR
            the direction, either left or right
        lineage : str
            the linege, which is the name of the pendant node
        time : float
            the time stamp of the event
        """
        if isinstance(chromosome, LinearChromosome):
            raise NotImplementedError('Inversion within intergenes not implemented'
                                      ' for linear chromosomes.')
        assert isinstance(chromosome, CircularChromosome)

        try:
            r = chromosome.get_affected_region(c1, c2, d)
        except CoordinateChoiceError:
            return False

        gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r

            # Get lengths from last intergene before modifying chromosome.
        specificlen = chromosome.intergenes[-1].specific_flanking[1]
        totallen = chromosome.intergenes[-1].total_flanking[1]

        if d == T_DIR.LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        sleftlen, srightlen, tleftlen, trightlen = 0, 0, 0, 0
        if c1 > c2:                     #The inversion wraps:
            sleftlen, srightlen, tleftlen, trightlen = \
                 chromosome.inversion_wrap_lengths(gpositions)

        inv = Inversion(int1, int2, c1, c2, specificlen, totallen, sleftlen,
                        tleftlen, srightlen, trightlen, lineage, time)

        #Nucleotide level inversion event
        chromosome.event_history.append(inv)

        #Gene-order level inversion event
        event = EventTwoBreakpoints(gpositions, chromosome, INV, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)

        segment = chromosome.obtain_segment(gpositions)
        chromosome.invert_segment(gpositions, igpositions)

        scar1 = chromosome.intergenes[igpositions[0]]
        scar2 = chromosome.intergenes[igpositions[-1]]

        scar1.length = leftlengths[0] + rightlengths[0]
        scar2.length = rightlengths[1] + leftlengths[1]

        assert scar1.length == len(inv.afterL)
        assert scar2.length == len(inv.afterR)      

        for gene in segment:
            self.all_gene_families[gene.family].append_event(time, INV, ";".join(map(str,[lineage, gene.gene_id])))

        return True


    def make_transposition(self, p: float, lineage: str, time: float):
        """
        Do a transposition on the given lineage. A segment is put in a new
        location chosen uniformly at random.
        """
        chromosome = self.node_genomes[lineage].select_random_chromosome()
        assert isinstance(chromosome, CircularChromosome)
        #Get a range of consecutive genome indices
        affected_indices = chromosome.obtain_affected_indices(p)
        #Get the corresponding segment
        segment = chromosome.obtain_segment(affected_indices)
        #Do the transposition
        old_pos, _ = chromosome.cut_and_paste(segment, affected_indices)

        event = EventThreeBreakpoints(affected_indices, old_pos, chromosome,
                                      POS, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)

        for gene in segment:
            self.all_gene_families[gene.family].append_event(time, POS, ";".join(map(str,[lineage, gene.gene_id])))


    def make_transposition_intergenic(self, chromosome: Chromosome,
                                      c1, c2, d: T_DIR, c3, lineage, time):
        if isinstance(chromosome, LinearChromosome):
            raise NotImplementedError('Transposition within intergenes not '
                                      'implemented for linear chromosomes.')

        try:
            r = chromosome.get_affected_region(c1, c2, d)
        except CoordinateChoiceError:
            return False

        gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r

        segment = chromosome.obtain_segment(gpositions)
        intergene_segment = chromosome.obtain_intergenic_segment(igpositions[1:])

        hereint, gene_position = chromosome.get_location_from_coord(c3, True)

        scar1 = chromosome.intergenes[igpositions[0]]    #stays put
        scar2 = chromosome.intergenes[hereint.position]  #left ig after transp
        scar3 = chromosome.intergenes[igpositions[-1]]   #right ig after transp
        assert scar3 != scar1 and scar2 != scar1, "segment can't be placed next to itself"

        new_segment = list()
        transposed_intergenes = list()

        # Get old lengths from last intergene before modifying chromosome.

        specificlen = chromosome.intergenes[-1].specific_flanking[1]
        totallen = chromosome.intergenes[-1].total_flanking[1]
        numintergenes = len(chromosome.intergenes)

        # If we insert in the intergene i, the gene must occupy the position i - 1
        # We store it for reference

        left_gene = chromosome.genes[hereint.position]

        # Now we pop the genes

        for gene in segment:
            new_segment.append(chromosome.genes.pop(chromosome.genes.index(gene)))

        # And now we insert the genes at the right of the gene we saved before

        position = chromosome.genes.index(left_gene) + 1

        for i, gene in enumerate(new_segment):
            chromosome.genes.insert(position + i, gene)

        # We move the intergene on the right also

        # We save the position for insertion

        here_intergene = chromosome.intergenes[hereint.position]

        # Remove the intergene segment (all intergenes except left breakpoint)

        for intergene in intergene_segment:
            transposed_intergenes.append(chromosome.intergenes.pop(chromosome.intergenes.index(intergene)))

        # And now we insert the transposed intergenes to the right of the insert point

        position = chromosome.intergenes.index(here_intergene) + 1

        for i, intergene in enumerate(transposed_intergenes):
            chromosome.intergenes.insert(position + i, intergene)

        # Finally, we modify the segments so that they have the right length

        herelengths = (c3 - hereint.sc1, hereint.sc2 - c3)

        if d == T_DIR.LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        if scar1 == scar3:    #Translocated segment placed in right intergene
                              #(this will never get called unless we allow such a thing)
            scar2.length = leftlengths[1] + (herelengths[0] - rightlengths[0]) + leftlengths[0]
            scar3.length = rightlengths[0] + herelengths[1]
            scar1 = scar2
        elif scar1 == scar2:  #Translocated segment placed in left intergene
                              #(this will never get called unless we allow such a thing)
            scar2.length = leftlengths[1] + herelengths[0]
            scar3.length = rightlengths[0] + (herelengths[1] - leftlengths[1]) + rightlengths[1]
            scar1 = scar3
        else:                 #Translocated segment placed in other intergene
            scar1.length = leftlengths[0] + rightlengths[1]
            scar2.length = leftlengths[1] + herelengths[0]
            scar3.length = rightlengths[0] + herelengths[1]

        trans = Transposition(int1, int2, c1, c2, hereint, c3, numintergenes,
                              specificlen, totallen, lineage, time)

        #Chromosome level origination event
        chromosome.event_history.append(trans)

        #Gene-order level origination event
        event = EventThreeBreakpoints(gpositions, gene_position, chromosome,
                                      POS, lineage, time)
        #chromosome.geneorder_history.append(event)
        self.geneorder_events[lineage].append(event)

        assert len(trans.afterH) == scar1.length
        assert len(trans.afterL) == scar2.length
        assert len(trans.afterR) == scar3.length

        for i, gene in enumerate(segment):
            self.all_gene_families[gene.family].append_event(time, POS, ";".join(map(str,[lineage, gene.gene_id])))

        return True


    def make_rewiring_edge(self, lineage, time):

        chromosome = self.node_genomes[lineage].select_random_chromosome()
        assert isinstance(chromosome, CircularChromosome)
        position = chromosome.select_random_position()
        interactome = self.node_genomes[lineage].interactome
        normalized_weights = af.normalize([d + 1 for n, d in interactome.degree()])  #type: ignore
        n1 = chromosome.genes[position]
        n2 = G_NPRNG().choice(sorted(interactome.nodes), 1, p=normalized_weights)[0]   #type: ignore

        while (str(n1) == str(n2)):
            n2 = G_NPRNG().choice(sorted(interactome.nodes), 1, p=normalized_weights)[0]  #type: ignore

        self.node_genomes[lineage].interactome.add_edge(str(n1), str(n2))
        self.all_gene_families[n1.family].append_event(time, "RW", ";".join(map(str, [lineage, n1, n2])))

    def make_remove_edge(self, lineage, time):

        myedges = list(self.node_genomes[lineage].interactome.edges())

        if len(myedges) == 0:
            return None

        myedge = myedges[G_RNG().randint(0, len(myedges) - 1)]
        self.node_genomes[lineage].interactome.remove_edge(*myedge)

        self.all_gene_families[myedge[0].split("_")[0]].append_event(time, "RM", ";".join([lineage, myedge[0], myedge[1]]))


    #def get_gene_family_tree(self):

    #    if len(self.gene_family["Gene_tree"].get_leaves()) < 3:
    #        return "None"
    #    else:
    #        return self.gene_family["Gene_tree"].write(format=1)


    def select_advanced_length(self, lineage: str, p: float, reps=100) \
        -> tuple[CircularChromosome, int, int, T_DIR]:
        """
        Return a pair of specific coordinates for intergenic regions according
        to `p` on a chromosome of `lineage`. The event must cover at least one gene

        Parameters
        ----------
        lineage: str
            the lineage to choose the chromosome from (the pendant node in 
            the species tree, e.g. n34)
        p: float
            1/p should be the expected (in nucleotides) difference between
            sc1 and sc2. in other words the expected number of intergenic
            nucleotides between the two breakpoints.
        reps: int
            try this many times to get a legal breakpoint pair.

        Returns
        -------
        Tuple[CircularChromosome, int, int, T_DIR]
            (sc1, sc2, direction) where sc1 and sc2 are specific intergenic
            breakpoint coordinates, meant to be breakpoints, and direction is
            one of {LEFT, RIGHT} indicating if sc2 is left or right of sc1.
        """
        chromosome = self.node_genomes[lineage].select_random_chromosome()
        assert isinstance(chromosome, CircularChromosome)

            #The total number of intergenic nucleotides can be retrieved from
            #the last intergenic location:
        assert chromosome.map_of_locations[-1].isIntergenic()
        intergenic_specific_length = chromosome.map_of_locations[-1].sc2

        success = False
        counter = 0
        while counter <= reps and success == False:
            counter += 1

            sc1 = chromosome.select_random_coordinate_in_intergenic_regions()
            d = G_NPRNG().choice((T_DIR.LEFT, T_DIR.RIGHT), p=[0.5, 0.5]) #type: ignore

            extension = G_NPRNG().geometric(p)

            if d == T_DIR.RIGHT:

                if sc1 + extension > intergenic_specific_length:
                    sc2 = sc1 + extension - intergenic_specific_length - 1
                    if sc2 < sc1:       # The event wraps to the right and
                        success = True  # doesn't cover the whole genome
                else:
                    sc2 = sc1 + extension
                    success = True

            else:
                assert d == T_DIR.LEFT

                if sc1 - extension < 0:
                    sc2 = intergenic_specific_length - (extension - sc1) + 1
                    if sc1 < sc2:       # The event wraps to the left and
                        success = True  # doesn't cover the whole genome
                else:
                    sc2 = sc1 - extension
                    success = True

            if success:

                l1 = chromosome.get_location_from_coord(sc1, True)
                l2 = chromosome.get_location_from_coord(sc2, True)
                if l1 != l2:
                    return chromosome, sc1, sc2, d

                # Verify that the event covers at least one gene

                r = chromosome.get_affected_region(sc1, sc2, d)
                
                if r == None:
                    success = False
                else:

                    genepositions, intergenepositions, leftlengths, rightlengths, int1, int2 = r
                            
                    if l1 != l2 and len(genepositions) != 0:
                        return chromosome, sc1, sc2, d
                    else:
                        success = False

        raise(CoordinateChoiceError)


    def return_cuts_by_event(self, event):
        """
        Return all the cuts of the event
        """
        if event.etype == ORIG:
            cuts = [event.sbp]
        elif event.etype == POS:
            cuts = [event.sbpL, event.sbpR, event.sbpH]
        elif event.etype == "N": # Transfer donor
            cuts = [event.sbpL, event.sbpR]
        elif event.etype == FER: # Transfer recipient
            cuts = [event.receptorsbp]
        else:
            cuts = [event.sbpL, event.sbpR]

        return cuts


    def propagate_cut(self, cut: int, event: GenomeCoordEvent) \
        -> tuple[bool, int, GenomeCoordEvent|None]:

        current_lineage: TreeNode = self.complete_tree & event.lineage   #ETE3's (bizarre) syntax for finding a node in a subtree
        chromosome = [chromosome for chromosome in self.node_genomes[current_lineage.name]][0]
        reversed_event_history = list(reversed(chromosome.event_history))
        index = reversed_event_history.index(event)
        
        finished = False
        adjust_index = True

        event2 = None
        
        while finished == False:

            if current_lineage == self.complete_tree: # We are at the root
                finished = True # We propagate the last time and then we end

            if adjust_index:
                reversed_event_history = reversed_event_history[index+1:]
                adjust_index = False

            for event2 in reversed_event_history:
                
                if event2.etype == "N": # Donor does not change the coordinates
                    pass

                elif event2.etype == FER:

                    assert isinstance(event2, Transfer)
                    lineage, cut = event2.afterToBeforeS_lineage(cut)

                    if lineage == current_lineage.name:
                        # I can remain in the same branch  
                        pass

                    else: # We need to propagate through a different branch

                        current_lineage = self.complete_tree & lineage
                        chromosome = [chromosome for chromosome in self.node_genomes[current_lineage.name]][0]
                        reversed_event_history = list(reversed(chromosome.event_history))
                        index = reversed_event_history.index(event2.sister_event)
                        adjust_index = True
                        break

                elif(event2.etype == LOSS and
                     (isinstance(event2, Loss) and event2.pseudogenize == True)):

                    # We are in a pseudogenized region.
                    # Four things can happen. 
                    # 1. The cut falls outside the pseudogenized region. We just continue passing the cut
                    # 2. The cut falls inside the pseudogenized region, in a previous breakpoint
                    # 3. The cut falls in the middle of a gene
                    # 4. The cut falls in the middle of an intergene (no divisions)

                    try:
                        cut = event2.afterToBeforeS(cut) # Case 1
    
                    except MapPseudogeneError:
                        
                        assert isinstance(event2, Loss)
                        piece, _ = event2.returnPieceAndCut(cut)
                        if type(piece) == Intergene:
                            # The cut has fallen into an intergene, we can keep propagating
                            pass
                        
                        if type(piece) == Gene:
                            #print("Now the cut is %s" % cut)
                            return False, cut, event2             
                        else:
                            pass
                            #print("Case 2")
                            #rint("The cut falls into a previous breakpoint, I keep going")          

                else:
                    cut = event2.afterToBeforeS(cut)
                
                    
             
            if finished == False and adjust_index == False: # Adjust index is true only when there has been a change
                                                            # to a different branch through a transfer event
            
                assert current_lineage.up
                current_lineage = current_lineage.up 
                chromosome = [chromosome for chromosome in self.node_genomes[current_lineage.name]][0]
                reversed_event_history = list(reversed(chromosome.event_history))

        #assert event2 is not None
        return True, cut, event2


    def return_all_events(self) -> list[GenomeCoordEvent]:
        
        all_events = list()
        
        for node in self.complete_tree.traverse("postorder"):   #type: ignore
                
            genome = self.node_genomes[node.name]                        
            chromosome = genome.chromosomes[0]
            chromosome.update_flankings()
            
            all_events += chromosome.event_history
        
        all_events = list(reversed(sorted(all_events, key=lambda x: x.time)))

        return all_events
        

    def init_divisions(self):
        """
        Initialize the divisions at the root
        """
        # First, we create a list with all the events ordered by time
        
        all_events = self.return_all_events()
        
        ########

        # Second, we traverse the events until the beginning

        initial_cuts: set[int] = set()
        
        self.gene2pseudogenecuts = dict()
        
        for event1 in all_events:
            cuts = self.return_cuts_by_event(event1)
            
            for cut in cuts:
                until_the_beginning, propagated_cut, p_event = self.propagate_cut(cut, event1)
                #print("We start with cut %s and propagate to cut %s" % (cut, propagated_cut))

                if until_the_beginning == True:
                    #print("Cut successfuly propagated", cut, "--->",propagated_cut, "Event:", event1.etype)
                    initial_cuts.add(propagated_cut)

                else:
                    # The cut has been propagated to a gene
                    assert isinstance(p_event, Loss)
                    gene, cut_within_gene = p_event.returnPieceAndCut(propagated_cut)

                    print(f"THIS SHOULDN'T WORK?! the type of gene is {type(gene)}")
                    gene_name = gene.family + "_" + gene.species + "_" + str(gene.gene_id) #TODO: the type of gene is Intergene?
                    if gene_name not in self.gene2pseudogenecuts:
                        self.gene2pseudogenecuts[gene_name] = set()
                    self.gene2pseudogenecuts[gene_name].add(cut_within_gene)

                    #print("Cut successfuly propagated within Gene", cut, "--->",cut_within_gene, "Event:", event1.etype, "Gene", gene_name)


        initial_chromosome = self.initial_genome.chromosomes[0]   
        all_cuts: set[int] = set()
        self.natural_cuts: list[tuple[int, int]] = list()

        # These are the natural cuts from the intergenes (the limits with the
        # genes):
        for intergene in initial_chromosome.iter_intergenes():

            cut1, cut2 = intergene.specific_flanking
            self.natural_cuts.append((cut1, cut2))
            
            all_cuts.add(cut1)
            all_cuts.add(cut2)


        all_cuts |= initial_cuts

        # These are the cuts surrounding the Genes. We don't need to treat these
        cuts_to_ignore = {(x1[1], x2[0]) for x1, x2 in zip(self.natural_cuts,
                                                           self.natural_cuts[1:] +
                                                           [self.natural_cuts[0]])}


        sorted_cuts =  sorted(all_cuts)
        initial_specific_flankings = zip(sorted_cuts,
                                         sorted_cuts[1:] + [sorted_cuts[0]])
        self.division_fam_id = 0

        self.initial_divisions: list[tuple[int, int]] = list()  # For debugging purposes

        for initial_specific_flanking in initial_specific_flankings:
            # initial_specific_flanking is a tuple (c1, c2)
            # We need to ignore the cuts where c1 is the right most extreme of
            # an intergene and c2 is the left most extreme of the next intergene
            if initial_specific_flanking in cuts_to_ignore:
                continue

            self.initial_divisions.append(initial_specific_flanking)
            #("Initial division", initial_specific_flanking)

            self.division_fam_id += 1

            intergene = initial_chromosome.return_intergene_by_coordinate(initial_specific_flanking[0])
            intergene.create_division(1, self.division_fam_id, initial_specific_flanking) # Identifier is 1 in the beginning
            division_family = DivisionFamily(self.division_fam_id, initial_specific_flanking)
            division_family.append_event(0, ORIG, "Root") # We register the origination

            self.all_division_families[self.division_fam_id] = division_family


    def redo_events_for_divisions(self) -> dict[str, Genome]:
        """
        Rerun the events on the tree, creating the pieces accordingly.
        """
        self.node_genomes_pieces: dict[str, Genome] = {} #:Pieces for genomes at each node
        all_genomes_pieces: dict[str, Genome] = {}       #:Pieces for all genomes
        lineage_counter = Counter()         #Number of active lineages per species

        self.gene_families_second = self.initial_gene_families
        self.initial_genome.init_pieces()
        self.node_genomes_pieces["Initial"] = self.initial_genome
        self.node_genomes_pieces["Root"] = copy.deepcopy(self.initial_genome)
        if self.save_all:
            all_genomes_pieces["Initial-0"] = self.initial_genome

        # We create a list of all the events (T and G) that we will order by time

        all_events: list[tuple[str, tuple[float, str|GenomeCoordEvent, str|Chromosome]]] = \
            [("T", x) for x in self.tree_events]

        for node in self.complete_tree.traverse():          #type: ignore
           genome = self.node_genomes[node.name]           
           for chromosome in genome:
               for event in chromosome.event_history:              
                   if event.etype == "N":
                       continue
                   assert isinstance(event, GenomeCoordEvent)
                   all_events.append(("G", (event.time, event, chromosome)))            

        # We start with the initial genome, that we preserved in a copy of the initial genome in the main simulation

        for items in sorted(all_events, key=lambda x: x[1][0]): # This sorts the events by the time
            
            # We unpack the events

            if items[0] == "T": # Tree event
                time, etype, lineages = items[1] # In the case that it is a species level event
                assert(isinstance(lineages, str))

                # Species level events
 
                if etype == "S":
                    pn, c1, c2 = lineages.split(";") # Parent, child1, child2
                    if self.save_all:
                        self.add_pieces(pn, all_genomes_pieces,
                                        lineage_counter, True)

                    self.make_speciation_divisions(time, pn, c1, c2)

                if etype == "E":
                    if self.save_all:
                        self.add_pieces(lineages, all_genomes_pieces,
                                        lineage_counter, True)

                    self.make_extinction_divisions(time, lineages)

                if etype == "F":               
                    if self.save_all:
                        self.add_pieces(lineages, all_genomes_pieces,
                                        lineage_counter, True)

                    self.make_end_divisions(time, lineages)
 

            else:               # Genome event
                time, event, _ = items[1]
                assert isinstance(event, GenomeCoordEvent)
                etype = event.etype
            
                # Genome level events

                if etype == TDUP:
                    assert isinstance(event, TandemDup)
                    lineage = self.make_tandemdup_divisions(time, event)
                elif etype == FER:
                    lineage, recipient = self.make_transfer_divisions(time, event)
                    if self.save_all:
                        self.add_pieces(recipient, all_genomes_pieces,
                                        lineage_counter)
                elif etype == LOSS:
                    assert isinstance(event, Loss)
                    lineage = self.make_loss_divisions(time, event)
                elif etype == INV:
                    assert isinstance(event, Inversion)
                    lineage = self.make_inversion_divisions(time, event)
                elif etype == POS:
                    assert isinstance(event, Transposition)
                    lineage = self.make_transposition_divisions(time, event)
                elif etype == ORIG:
                    assert isinstance(event, Origination)
                    lineage = self.make_origination_divisions(time, event)
                else:
                    raise NotImplementedError(f"Unkown event type {etype}.")

                if self.save_all:
                    self.add_pieces(lineage, all_genomes_pieces, lineage_counter)
    
        return all_genomes_pieces


    def make_speciation_divisions(self, time, pn, c1, c2):
        
        genome_pn: Genome = self.node_genomes_pieces[pn]
        genome1 = Genome()
        genome2 = Genome()

        self.node_genomes_pieces[c1] = genome1
        self.node_genomes_pieces[c2] = genome2

        for chromosome in genome_pn:

            ch1 = CircularChromosome(chromosome.name)
            ch2 = CircularChromosome(chromosome.name)

            genome1.chromosomes.append(ch1)
            genome2.chromosomes.append(ch2)
 
            for piece in chromosome.pieces: # We iterate the pieces in the parent chromosome

                if isinstance(piece, Gene):

                    # We need to insert new genes in the children chromosomes. We need to keep track of the length, the orientation and the total coordinates

                    gene = piece # For the sake of clarity

                    new_id1 = self.return_new_gene_ids_for_segment_with_divisions([gene])[0]
                    new_id2 = self.return_new_gene_ids_for_segment_with_divisions([gene])[0]

                    new_gene1 = Gene()
                    new_gene2 = Gene()

                    new_gene1.gene_id = new_id1
                    new_gene2.gene_id = new_id2

                    new_gene1.ptype = gene.ptype
                    new_gene2.ptype = gene.ptype

                    new_gene1.orientation = gene.orientation
                    new_gene2.orientation = gene.orientation
 
                    new_gene1.length = gene.length
                    new_gene2.length = gene.length

                    new_gene1.total_flanking = gene.total_flanking
                    new_gene2.total_flanking = gene.total_flanking

                    new_gene1.family = gene.family
                    new_gene2.family = gene.family

                    new_gene1.species = c1
                    new_gene2.species = c2

                    ch1.pieces.append(new_gene1)
                    ch2.pieces.append(new_gene2)

                if isinstance(piece, Division):

                    division = piece

                    # If the piece is a division, I need to keep track also of the identity within the gene family
                    new_identity1 = self.all_division_families[division.family].obtain_new_identifier()
                    new_identity2 = self.all_division_families[division.family].obtain_new_identifier()

                    division1 = Division(new_identity1, division.family)
                    division2 = Division(new_identity2, division.family)

                    division1.total_flanking = division.total_flanking
                    division2.total_flanking = division.total_flanking
                    
                    division1.specific_flanking = division.specific_flanking
                    division2.specific_flanking = division.specific_flanking
                    
                    division1.orientation = division.orientation
                    division2.orientation = division.orientation
                    
                    division1.length = division.length
                    division2.length = division.length
                    
                    division1.ptype = division.ptype
                    division2.ptype = division.ptype

                    division1.species = c1
                    division2.species = c2
                    
                    nodes = [pn,
                             division.identity,
                             c1,
                             division1.identity,
                             c2,
                             division2.identity
                            ]
                    self.all_division_families[division.family].append_event(time, "S", ";".join(map(str, nodes)))
                    
                    ch1.pieces.append(division1)
                    ch2.pieces.append(division2)
        
        ch1.update_coordinates()   
        ch2.update_coordinates()   
 

    def make_extinction_divisions(self, time, lineage):
        
        chromosome = [x for x in self.node_genomes_pieces[lineage]][0] 

        for piece in chromosome.pieces:
            if isinstance(piece, Division):
                self.all_division_families[piece.family].append_event(time, "E", ";".join(map(str,[lineage, piece.identity])))


    def make_end_divisions(self, time, lineage):
        
        chromosome = [x for x in self.node_genomes_pieces[lineage]][0]

        for piece in chromosome.pieces:
            if isinstance(piece, Division):
                self.all_division_families[piece.family].append_event(time, "F", str(lineage) + ";" + str(piece.identity)) 
    

    def select_pieces(self, chromosome: Chromosome, tcL, tcR) \
        -> tuple[list[Gene | Division], list[int], bool]:
        """
        Select the pieces affected by an event given the total coordinates of
        the breakpoints.

        Returns
        -------
        tuple[list[Gene | Division], list[int], bool]
            (pieces_affected, indices_affected, wrapping) where pieces_affected
            is the list of pieces affected (genes and divisions),
            indices_affected is the list of indices of the pieces affected, and
            wrapping is a boolean indicating whether the event wraps around the
            circular chromosome.
        """
       
        start = False
        end = False

        pieces_affected: list[Gene | Division] = []
        indices_affected = []

        # If the event affects at the end or the beginning of the chromosome:

        wrapping = False

        if tcL == chromosome.pieces[-1].total_flanking[1]:
            tcL = 0
            wrapping = True

        if tcR == 0:
            tcR = chromosome.pieces[-1].total_flanking[1]
            wrapping = True

        # Cycle through the pieces twice in case the event wraps around
        for index, piece in enumerate(itertools.cycle(chromosome.pieces)):               

            pfL, pfR = piece.total_flanking

            if pfL == tcL:
                start = True
            if start == True:
                if index >= len(chromosome.pieces): # We are in the second cycle
                    indices_affected.append(index - len(chromosome.pieces))
                    wrapping = True
                else:
                    indices_affected.append(index)
                pieces_affected.append(piece)
            if pfR == tcR and start == True:
                end = True
            if end == True:
                break

            if index >= 2 * len(chromosome.pieces) + 1: # FIX This is here just for debugging purposes
                raise(Exception(f'Piece index cannot be found "{tcL, tcR}"'))

        return pieces_affected, indices_affected, wrapping


    def make_tandemdup_divisions(self, time: float, event: TandemDup):
        
        lineage = event.lineage        
        chromosome = [x for x in self.node_genomes_pieces[lineage]][0] 
        
        tcL = event.tbpL
        tcR = event.tbpR

        pieces_to_duplicate, indices_to_duplicate, _ = self.select_pieces(chromosome, tcL, tcR)

        # We copy the pieces

        pieces_duplicated = list()

        # We update the identifiers of the duplicated and the not duplicated pieces

        insert_index = indices_to_duplicate[-1] + 1

        # The gene identifiers need to be assigned in this ordered
        # to make it coincide with the forward simulation

        new_gene_identifiers1 = self.return_new_gene_ids_for_segment_with_divisions(pieces_to_duplicate)
        new_gene_identifiers2 = self.return_new_gene_ids_for_segment_with_divisions(pieces_to_duplicate)
        
        #####

        i = 0 # This is to keep track of the genes

        for original_piece in pieces_to_duplicate:

            duplicated_piece = copy.deepcopy(original_piece)
            
            if isinstance(original_piece, Division):

                division_family = original_piece.family
                parent_id = original_piece.identity
                new_id1 = self.all_division_families[division_family].obtain_new_identifier()
                new_id2 = self.all_division_families[division_family].obtain_new_identifier()            
                
                original_piece.identity = new_id1
                assert isinstance(duplicated_piece, Division)
                duplicated_piece.identity = new_id2        

                self.all_division_families[division_family].append_event(time, TDUP, ";".join(map(str,[lineage, parent_id, lineage, new_id1, lineage, new_id2])))
            
            else:
                assert isinstance(original_piece, Gene) and isinstance(duplicated_piece, Gene)
                original_piece.gene_id = new_gene_identifiers1[i]
                duplicated_piece.gene_id = new_gene_identifiers2[i]
                i+=1                
            
            pieces_duplicated.append(duplicated_piece)


        chromosome.pieces = chromosome.pieces[0:insert_index] + pieces_duplicated + chromosome.pieces[insert_index:] 
        chromosome.update_coordinates()

        return lineage


    def make_transfer_divisions(self, time, event) -> tuple[str, str]:
        """
        Do a horizontal transfer event, while updating intergenic regions.

        Returns
        -------
        tuple[str, str]
            (donor_lineage, recipient_lineage)
        """
        donor_lineage = event.donorlineage
        recipient_lineage = event.receptorlineage        

        donor_chromosome = [x for x in self.node_genomes_pieces[donor_lineage]][0] 
        recipient_chromosome = [x for x in self.node_genomes_pieces[recipient_lineage]][0] 
        
        tcL = event.tbpL
        tcR = event.tbpR
        
        pieces_to_transfer, _, _ = self.select_pieces(donor_chromosome, tcL, tcR)

        insert_after_this_piece = None                

        insertion_point = event.receptortbp

        for piece in itertools.cycle(recipient_chromosome.pieces):
            _, pfR = piece.total_flanking
            #print(pfL, pfR, insertion_point)
            
            if pfR == insertion_point: 
                insert_after_this_piece = piece
           
                break

        

        # The gene identifiers need to be assigned in this ordered
        # to make it coincide with the forward simulation

        new_gene_identifiers1 = self.return_new_gene_ids_for_segment_with_divisions(pieces_to_transfer)
        new_gene_identifiers2 = self.return_new_gene_ids_for_segment_with_divisions(pieces_to_transfer)

        i = 0 # This is to keep track of the genes
        pieces_transferred = list()

        for original_piece in pieces_to_transfer:

            transferred_piece = copy.deepcopy(original_piece)
            
            if isinstance(original_piece, Division):

                division_family = original_piece.family
                parent_id = original_piece.identity
                new_id1 = self.all_division_families[division_family].obtain_new_identifier()
                new_id2 = self.all_division_families[division_family].obtain_new_identifier()            
                original_piece.identity = new_id1
                assert isinstance(transferred_piece, Division)
                transferred_piece.identity = new_id2     
                           
                self.all_division_families[division_family].append_event(time, FER, ";".join(map(str,[donor_lineage, parent_id, donor_lineage, new_id1, recipient_lineage, new_id2]))) 
                
            
            else:
                assert isinstance(original_piece, Gene) and isinstance(transferred_piece, Gene)
                original_piece.gene_id = new_gene_identifiers1[i]
                transferred_piece.gene_id = new_gene_identifiers2[i]
                transferred_piece.species = recipient_lineage
                i+=1                
            
            pieces_transferred.append(transferred_piece)


        assert insert_after_this_piece is not None
        insert_index = recipient_chromosome.pieces.index(insert_after_this_piece) + 1


        recipient_chromosome.pieces = recipient_chromosome.pieces[0:insert_index] + pieces_transferred + recipient_chromosome.pieces[insert_index:]

        #recipient_chromosome.print_pieces()
        recipient_chromosome.update_coordinates()

        return donor_lineage, recipient_lineage


    def make_loss_divisions(self, time: float, event: Loss):

        lineage = event.lineage        
        chromosome = [x for x in self.node_genomes_pieces[lineage]][0] 
        pseudo = event.pseudogenize

        
        tcL = event.tbpL
        tcR = event.tbpR
        
        if pseudo and event.wraps(): # These are the right coordinates if there is a pseudogenization
                                     # and the event wraps
            tcL = event.tc1
            tcR = event.tc2

        pieces_to_lose, indices_to_lose, wrapping = self.select_pieces(chromosome, tcL, tcR)

        if not pseudo:
            
            chromosome.pieces = [piece for piece in chromosome.pieces if piece not in pieces_to_lose]
            for piece in pieces_to_lose:
                if isinstance(piece, Division):
                    self.all_division_families[piece.family].append_event(time, LOSS, ";".join(map(str,[lineage, piece.identity])))
        else:

            replacements = dict()

            for index, piece in zip(indices_to_lose, pieces_to_lose):
                if isinstance(piece, Gene):

                    replacements[piece] = list()

                    gene_name = piece.family + "_" + piece.species + "_" + str(piece.gene_id)

                    cuts = {0, piece.length}

                    #print("The cuts in the gene are", self.gene2pseudogenecuts)
                    #print(gene_name)
                    
                    if gene_name in self.gene2pseudogenecuts:
                        cuts = cuts.union(self.gene2pseudogenecuts[gene_name])
                    cuts = sorted(list(cuts))

                    # We need to make as many divisions as cuts + 1
                    

                    for cut1, cut2 in zip(cuts, cuts[1:]):
                    
                        # Insert new family

                        self.division_fam_id += 1

                        division_family = DivisionFamily(self.division_fam_id, (0,0)) 
                        division_family.initial_orientation = piece.orientation
                        division_family.append_event(time, ORIG, lineage) # We register the origination. 
                        division = Division(1, self.division_fam_id, (0,0))
                        division.length = cut2 - cut1
                        division.total_flanking = (piece.total_flanking[0] + cut1, piece.total_flanking[0] + cut2)
                        division.species = lineage
                        replacements[piece].append(division)
                        division.initial_sequence = gene_name

                        self.all_division_families[self.division_fam_id] = division_family

            for gene, divisions in replacements.items():
                insert = (chromosome.pieces).index(gene)
                (chromosome.pieces).pop(insert)
                for i, division in enumerate(divisions):
                    (chromosome.pieces).insert(insert + i, division)

                        
        if wrapping == True:
            while chromosome.pieces[0].ptype == "Divi":
                chromosome.pieces = chromosome.pieces[1:] + [chromosome.pieces[0]]

        # Now we start adjusting the pieces until there are no divisions in the beginning
        # and the unaffected pieces remain in the same index as before
            
        while (chromosome.pieces[0].ptype == "Divi"):
                piece = (chromosome.pieces).pop(0)
                (chromosome.pieces).append(piece)
            
        chromosome.update_coordinates()

        #chromosome.print_pieces()
        return lineage


    def make_transposition_divisions(self, time, event: Transposition):

        lineage = event.lineage        
        chromosome = [x for x in self.node_genomes_pieces[lineage]][0] 
        
        tcL = event.tbpL
        tcR = event.tbpR

        pieces_to_transpose, _, _ = self.select_pieces(chromosome, tcL, tcR)

        #genes = [piece for piece in chromosome.pieces if piece.ptype == "Gene"]

        # The tranposed pieces will be in position tbpH

        tbpH = event.tbpH
        insert_after_this_piece = None        


        for piece in itertools.cycle(chromosome.pieces):
            pfL, pfR = piece.total_flanking
            if pfR == tbpH: 
                insert_after_this_piece = piece
                break
        
        # We remove first the pieces to transpose

        chromosome.pieces = [piece for piece in chromosome.pieces
                             if piece not in pieces_to_transpose]

        # We insert the pieces

        assert insert_after_this_piece is not None
        insert_index = chromosome.pieces.index(insert_after_this_piece) + 1 
        chromosome.pieces = chromosome.pieces[0:insert_index] + pieces_to_transpose + chromosome.pieces[insert_index:] 
    
        # Now we start adjusting the pieces until there are no divisions in the beginning
        # and the unaffected pieces remain in the same index as before
            
        while (chromosome.pieces[0].ptype == "Divi"):
                piece = (chromosome.pieces).pop(0)
                (chromosome.pieces).append(piece)
        

        chromosome.update_coordinates()

        return lineage


    def make_inversion_divisions(self, time, event: Inversion):
        
        lineage = event.lineage
        chromosome = next(iter(self.node_genomes_pieces[lineage]))
        
        tcL = event.tbpL
        tcR = event.tbpR
        
        pieces_to_invert, indexes_to_invert, wrapping = self.select_pieces(chromosome, tcL, tcR)

        genes = [piece for piece in chromosome.pieces if isinstance(piece, Gene)]
        gene2index = {gene:index for index, gene in enumerate(genes)}

        pieces_to_invert = list(reversed(copy.deepcopy(pieces_to_invert)))

        for index, replacement in zip(indexes_to_invert, pieces_to_invert):
            chromosome.pieces[index] = replacement
            replacement.change_orientation()

        # Now we adjust the indexes if there has been a wrapping event

        if wrapping == True:

            # We search the index position of a gene not affected by the event

            gene_ref_index = 0
            gene_ref = None

            for piece in chromosome.pieces:
                if isinstance(piece, Gene):
                    gene_ref_index +=1
                    gene_ref = piece
                    if piece not in pieces_to_invert:
                        break

            # Now we start adjusting the pieces until there are no divisions in the beginning
            # and the unaffected pieces remain in the same index as before
            
            assert gene_ref is not None
            while ((chromosome.pieces[0].ptype == "Divi") or
                   chromosome.get_index_gene(gene_ref) != gene2index[gene_ref]):
                piece = (chromosome.pieces).pop(0)
                (chromosome.pieces).append(piece)

        chromosome.update_coordinates()

        return lineage


    def make_origination_divisions(self, time, event: Origination):

        lineage = event.lineage        
        chromosome = [x for x in self.node_genomes_pieces[lineage]][0] 
        
        tc = event.interval.specificToTotal(event.sbp)

        gene = Gene()
        gene.species = lineage
        gene.length = event.genelen
        gene_family_id = event.gene_family
        gene.family = gene_family_id

        gene.orientation = self.all_gene_families[gene.family].initial_orientation 

        for piece in itertools.cycle(chromosome.pieces):
            _, pfR = piece.total_flanking
            if pfR == tc:
                insert_after_this_piece = piece
                break

        gene_family = GeneFamily(event.gene_family, time)
        gene_family.length = int(gene.length)
        gene_family.genes.append(gene)
        gene.gene_id = gene_family.obtain_new_gene_id()

        self.gene_families_second[gene_family_id] = gene_family
        self.gene_families_second[gene.family].append_event(time, ORIG, lineage)

        # We insert the pieces

        assert isinstance(insert_after_this_piece, (Gene, Division))
        insert_index = chromosome.pieces.index(insert_after_this_piece) + 1 
        chromosome.pieces.insert(insert_index, gene)
        chromosome.update_coordinates()

        return lineage


    def write_division_trees(self, division_tree_folder: Path):
        """
        This function writes the division trees and a file documenting the lenghts of every
        division in the initial genome
        """
        division_tree_folder.mkdir(parents=True, exist_ok=True)

        for division_family_name, division_family in self.all_division_families.items():
            complete_tree, pruned_tree, _ = division_family.generate_tree()

            with open(division_tree_folder / (str(division_family_name) + "_completetree.nwk"), "w") as f:
                f.write(complete_tree)

            if pruned_tree != None:
                with open(division_tree_folder / (str(division_family_name) + "_prunedtree.nwk"), "w") as f:
                    f.write(pruned_tree)

        with open(division_tree_folder / DIVISIONLENGTHS, "w") as f:
            for division_family_name, division_family in self.all_division_families.items():
                f.write("\t".join(list(map(str,[division_family_name, len(division_family)]))) + "\n")


def get_scaling_factor(scaled_file, events_file, scale) -> tuple[float, float]:
    """
    Compute the scaling factor given the scaled file and the events file.

    Returns
    -------
    tuple[float, float]
        (beginning_time, scaling_factor) where `beginning_time` is the time of
        the root speciation and `scaling_factor` is the factor to multiple the
        times in the events by.
    """
    # First I read where the root is:

    with open(scaled_file) as f:
        f.readline()
        for l in f:
            t, _, nodes = l.strip().split("\t")
            if float(t) == 0:
                eroot = nodes.split(";")[0]

    # Second, I read the total length of the tree
    with open(events_file) as f:                
        for l in f:                    
            t, _, nodes = l.strip().split("\t")
            node1 = nodes.split(";")[0]
            if node1 == eroot:
                beginning_time = float(t)
            continue
        t, _, nodes = l.strip().split("\t")

        totaltime = float(t) - beginning_time

    return beginning_time, scale / totaltime


def write_geneorder_events(outfile: Path, events: list[GeneOrderEvent],
                           scale=0.0, scaled_file: Path|None = None,
                           events_file: Path|None = None):
    header = "TIME\tEVENT\tBREAKPOINTS\tLENGTH\tCHROMOSOME\n"
    if scaled_file:
        beginning_time, mfactor = get_scaling_factor(scaled_file, events_file, scale)

    with open(outfile, "w") as f:
        f.write(header)
        for event in sorted(events, key=lambda x: x.time):
            etype, time, breakpoints, length, name = event.return_info()

            if scale: # Only working if Species Tree has been scaled to the same distance!
                time = (float(time) - beginning_time) * mfactor

            line = [str(time), etype, breakpoints, str(length), name]
            line = "\t".join(line) + "\n"
            f.write(line)



