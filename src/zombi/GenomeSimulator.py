from __future__ import annotations

import ete3
from ete3.coretype.tree import TreeNode
import numpy
import copy
import random
import os
import networkx as nx
from typing import Tuple, Set, Dict, Union

from . import AuxiliarFunctions as af
from .Events import Loss, Origination, TandemDup, Inversion, Transfer, Transposition
from .Genomes import Chromosome, CircularChromosome, Gene, GeneFamily, Genome, DivisionFamily, Intergene
from .Genomes import T_DIR, LEFT, RIGHT, Intergene, LinearChromosome

from Bio.SeqFeature import SeqFeature

class GenomeSimulator():
    """
    Attributes
    ----------
    active_genomes: Set[str]
        list of lineages that currently exist in the gene tree, each lineage is
        represented as a string indicating the pendant node (e.g. n37)
    all_genomes: Dict[str, Genome]
        map lineage name to Genome (gene order)
    """

    def __init__(self, parameters, events_file, root_genome: str=''):

        self.parameters = parameters

        mseed = self.parameters["SEED"]
        if mseed != 0:
            random.seed(parameters["SEED"])
            numpy.random.seed(parameters["SEED"])

        self.tree_events = self._read_events_file(events_file)
        self.distances_to_start = self._read_distances_to_start(events_file) # Only useful when computing assortative transfers
        self.complete_tree = self._read_tree(events_file.replace("Events.tsv", "CompleteTree.nwk"))

        self.all_genomes: Dict[str, Genome] = dict()
        self.all_gene_families: Dict[str, GeneFamily] = dict()

        self.all_division_families: Dict[str, DivisionFamily ] = dict()

        self.gene_families_counter = 0
        self.active_genomes: Set[str] = set()

        if self.parameters["RATE_FILE"] != "False":
            if self.parameters["SCALE_RATES"] == "True":
                self.crown_length = self._read_crown_length(events_file.replace("Events", "Lengths"))
                self.empirical_rates = af.read_empirical_rates(rates_file=self.parameters["RATE_FILE"], scale_rates= self.crown_length)
            else:
                self.empirical_rates = af.read_empirical_rates(rates_file=self.parameters["RATE_FILE"])

        self.root_genome_file = root_genome     #Get root genome from GFF file.
        if root_genome and not os.path.exists(root_genome):
            raise(Exception(f"Root genome file {root_genome} not found."))

        # A list to keep track of all the event coordinates

        self.event_coordinates = list()         #NOTE: remove this

    def write_genomes(self, genome_folder, intergenic_sequences = False):

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        for genome_name,genome in self.all_genomes.items():

            with open(os.path.join(genome_folder, genome_name + "_GENOME.tsv"), "w") as f:

                header = ["POSITION", "GENE_FAMILY", "ORIENTATION", "GENE_ID"]
                header = "\t".join(map(str, header)) + "\n"
                f.write(header)

                for chromosome in genome:
                    for index, gene in enumerate(chromosome):

                        line = [index, gene.gene_family, gene.orientation, gene.gene_id]
                        line = "\t".join(map(str,line)) +"\n"
                        f.write(line)

            if intergenic_sequences == True:

                with open(os.path.join(genome_folder, genome_name + "_LENGTHS.tsv"), "w") as f:

                    header = ["POSITION", "IDENTITY", "LENGTH"]
                    header = "\t".join(map(str, header)) + "\n"
                    f.write(header)

                    for chromosome in genome:
                        i = 0
                        for j, gene in enumerate(chromosome.genes):
                            line = [i, "G(" + str(gene.gene_family) + "_" + str(gene.gene_id) + ")", str(gene.length)]
                            line = "\t".join(map(str, line)) + "\n"
                            f.write(line)
                            i += 1
                            line = [i, "I", str(chromosome.intergenes[j].length)]
                            line = "\t".join(map(str, line)) + "\n"
                            f.write(line)
                            i += 1
        
    def write_division_coordinates(self, genome_folder):

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        for genome_name, genome in self.all_genomes_second.items():

            with open(os.path.join(genome_folder, genome_name + "_DIVISIONS.tsv"), "w") as f:

                header = ["DIVISION_FAMILY", "TYPE", "IDENTITY", "LENGTH", "FLANKING_LEFT_SPECIFIC", "FLANKING_RIGHT_SPECIFIC", "TOTAL_LEFT", "TOTAL_RIGHT", "ORIENTATION"]
                header = "\t".join(map(str, header)) + "\n"
                f.write(header)

                for chromosome in genome:
                    chromosome.obtain_locations()

                    for gene, intergene in zip(chromosome, chromosome.iter_intergenes()):

                        scl, scr = gene.total_flanking
                        tcl, tcr = gene.specific_flanking
                        line = "\t".join(list(map(str, [gene.gene_family, "GENE", gene.gene_id, gene.length, scl, scr, tcl, tcr, gene.orientation ]))) + "\n"
                        f.write(line)

                        for division in intergene:

                            scl = division.specific_flanking[0]
                            scr = division.specific_flanking[1]
                            tcl = chromosome.return_total_coordinate_from_specific_coordinate(scl)
                            tcr = chromosome.return_total_coordinate_from_specific_coordinate(scr)

                            line = "\t".join(list(map(str, [division.division_family, "DIVISION", division.identity, len(division), scl, scr, tcl, tcr, division.orientation ]))) + "\n"
                            f.write(line)


    def write_gene_family_info(self, genome_folder:str,
                               filename="GeneFamily_info.tsv"):
        """
        Write a TSV file with containing gene id, gff gene id, and start and
        end coordinates for every gene in `self.all_gene_families`.

        Parameters
        ----------
        genome_folder : str
            the folder
        filename : str, optional
            the filename to use, by default "GeneFamily_info.tsv"
        """
        with open(os.path.join(genome_folder, filename), "w") as f:
            header = ["GENE_FAMILY", "GFF_ID", "START", "END"]
            f.write("\t".join(map(str, header)) + "\n")

            for gene_family_name, gene_family in self.all_gene_families.items():
                if gene_family.gff_id:
                    line = "\t".join([gene_family_name, str(gene_family.gff_id),
                                      str(gene_family.genes[0].start+1),
                                      str(gene_family.genes[0].end)]) + "\n"
                    f.write(line)

    def write_gene_family_GFF_ids(self, genome_folder:str,
                               filename="GeneFamily_GFF_ids.tsv"):
        """
        Write a TSV file with containing gene id and gff gene id for every gene
        in `self.all_gene_families`.

        Parameters
        ----------
        genome_folder : str
            the folder
        filename : str, optional
            the filename to use, by default "GeneFamily_GFF_ids.tsv"
        """
        with open(os.path.join(genome_folder, filename), "w") as f:
            header = ["GENE_FAMILY", "GFF_ID"]
            f.write("\t".join(map(str, header)) + "\n")

            for gene_family_name, gene_family in self.all_gene_families.items():
                if gene_family.gff_id:
                    line = "\t".join([gene_family_name, str(gene_family.gff_id)]) + "\n"
                    f.write(line)
                
            
    def write_gene_family_lengths(self, genome_folder):

        with open(os.path.join(genome_folder, "GeneFamily_lengths.tsv"), "w") as f:
            header = ["GENE_FAMILY", "LENGTH"]
            header = "\t".join(map(str, header)) + "\n"
            f.write(header)

            for gene_family_name, gene_family in self.all_gene_families.items():
                line = "\t".join([gene_family_name, str(gene_family.length)]) + "\n"
                f.write(line)

    def write_gene_family_events(self, gene_family_events_folder):

        if not os.path.isdir(gene_family_events_folder):
            os.mkdir(gene_family_events_folder)

        for gene_family_name, gene_family in self.all_gene_families.items():

            with open(os.path.join(gene_family_events_folder, gene_family_name + "_events.tsv"),"w") as f:

                header = ["TIME","EVENT","NODES"]

                header = "\t".join(map(str, header)) + "\n"

                f.write(header)

                for time, event, nodes in gene_family.events:
                    
                    line = [time, event, nodes]

                    line = "\t".join(map(str, line)) + "\n"
                    f.write(line)
        

    def write_gene_trees(self, gene_tree_folder, gene_trees = True, reconciliations = False):

        if not os.path.isdir(gene_tree_folder):
            os.mkdir(gene_tree_folder)

        for gene_family_name, gene_family in self.all_gene_families.items():

            complete_tree, pruned_tree, rec_tree = gene_family.generate_tree()

            if gene_trees == True:

                with open(os.path.join(gene_tree_folder, gene_family_name + "_completetree.nwk"), "w") as f:
                    f.write(complete_tree)
                if pruned_tree != None:
                    with open(os.path.join(gene_tree_folder, gene_family_name + "_prunedtree.nwk"), "w") as f:
                        f.write(pruned_tree)
            if reconciliations == True:
                with open(os.path.join(gene_tree_folder, gene_family_name + "_rec.xml"), "w") as f:
                    f.write(rec_tree)

    def write_events_per_branch(self, events_per_branch_folder, scale, scaled_file, events_file):

        if not os.path.isdir(events_per_branch_folder):
            os.mkdir(events_per_branch_folder)

        events_per_branch = dict()

        for gene_family_name, gene_family in self.all_gene_families.items():

            for time, event, nodes in gene_family.events:

                name = nodes.split(";")[0]

                if name not in events_per_branch:
                    events_per_branch[name] = list()

                if event == "S" or event == "E" or event == "F":
                    continue

                elif event == "T":
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

                elif event == "D":

                    handle = nodes.split(";")
                    new_nodes = ";".join(map(lambda x: gene_family_name + "_" + x,[handle[1],handle[3],handle[5]]))
                    events_per_branch[name].append((time, event, new_nodes))

                else:

                    gene_id = nodes.split(";")[-1]
                    events_per_branch[name].append((time, event, gene_family_name + "_" + gene_id))


        for name, events in events_per_branch.items():

            with open(os.path.join(events_per_branch_folder, name + "_branchevents.tsv"), "w") as f:

                header = ["TIME", "EVENT", "NODES"]
                header = "\t".join(map(str, header)) + "\n"

                f.write(header)

                for time, event, nodes in sorted(events, key = lambda x: float(x[0])):

                    line = [str(time), event, nodes]
                    line = "\t".join(line) + "\n"
                    f.write(line)
        
            if scale != 0: # Only working if Species Tree has been scaled too the same distance!

                # First I read where the root is:

                with open(scaled_file) as f:
                    f.readline()
                    for l in f:
                        t, event, nodes = l.strip().split("\t")
                        if float(t) == 0:
                            eroot = nodes.split(";")[0]


                # Second, I read the total length of the tree

                with open(events_file) as f:                
                    for l in f:                    
                        t, event, nodes = l.strip().split("\t")
                        node1 = nodes.split(";")[0]
                        if node1 == eroot:
                            beginning_time = float(t)
                        continue
                    t, event, nodes = l.strip().split("\t")

                    totaltime = float(t) - beginning_time

                    mfactor = scale / totaltime


                with open(os.path.join(events_per_branch_folder, name + "_brancheventsscaled.tsv"), "w") as f:

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


        genome_names = [x for x in self.all_genomes.keys()]
        gene_family_names = [str(x) for x in self.all_gene_families.keys()]

        # For clarity, I start with Initial Genome

        genome_names[0], genome_names[1] = genome_names[1], genome_names[0]

        data = list()
        data.append(["GENOME"] + gene_family_names)


        for genome_name in genome_names:
            line = dict()

            genome = self.all_genomes[genome_name]

            for gene_family_name in gene_family_names:
                if gene_family_name not in line:
                    line[gene_family_name] = 0

            for chromosome in genome:
                for index, gene in enumerate(chromosome):
                    line[gene.gene_family] += 1

            data.append([genome_name] + [str(line[fm]) for fm in gene_family_names])

        # We will transpose the data

        mlenx = len(data)
        mleny = len(data[0])

        with open(os.path.join(profiles_folder, "Profiles.tsv"), "w") as f:
            for i in range(mleny):
                line = list()
                for j in range(mlenx):
                    line.append(str(data[j][i]))
                line = "\t".join(line) + "\n"
                f.write(line)


    def write_interactomes(self, genome_folder):

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        for genome_name, genome in self.all_genomes.items():

            if not hasattr(genome, "interactome"):
                continue

            with open(os.path.join(genome_folder, genome_name + "_INTERACTOME.tsv"), "w") as f:

                header = ["GENE_1", "GENE_2"]
                header = "\t".join(map(str, header)) + "\n"
                f.write(header)

                for g1, g2 in genome.interactome.edges:
                    f.write("\t".join([g1, g2]) + "\n")

    def write_family_rates(self, genome_folder):

        with open(os.path.join(genome_folder, "Family_rates.tsv"), "w") as f:
            header = ["GENE_FAMILY", "D", "T", "L"]
            header = "\t".join(map(str, header)) + "\n"
            f.write(header)

            for gene_family_name, gene_family in self.all_gene_families.items():

                d = gene_family.rates["DUPLICATION"]
                t = gene_family.rates["TRANSFER"]
                l = gene_family.rates["LOSS"]

                f.write("\t".join(map(str,[gene_family_name, d,t,l])) + "\n")


    def _read_events_file(self, events_file):

        events = list()
        with open(events_file) as f:
            f.readline()
            for line in f:
                handle = line.strip().split("\t")
                events.append(handle)
        return events

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

    def return_new_identifiers_for_segment(self, segment):

        new_identifiers = list()

        for gene in segment:
            gf = gene.gene_family
            new_id = self.all_gene_families[gf].obtain_new_gene_id()
            new_identifiers.append(new_id)

        return new_identifiers

    def return_new_identifiers_for_segment_with_divisions(self, segment):

        new_identifiers = list()

        for gene in segment:
            gf = gene.gene_family
            new_id = self.gene_families_second[gf].obtain_new_gene_id()
            new_identifiers.append(new_id)

        return new_identifiers

    def fill_genome(self, intergenic_sequences = False, family_rates = False, interactome = False):

        genome = Genome()
        genome.species = "Root"
        time = 0

        initial_genome_size = self.parameters["INITIAL_GENOME_SIZE"].split(";")
        shape = "C"

        for n_genes in initial_genome_size:

            if shape == "L":
                chromosome = LinearChromosome()
                chromosome.shape = "L"
            elif shape == "C":
                chromosome = CircularChromosome()
                chromosome.shape = "C"

            if intergenic_sequences == True:

                chromosome.has_intergenes = True
                mean_length = int(self.parameters["INTERGENE_LENGTH"])
                intergene_lengths = [int(x * mean_length * int(n_genes)) for x in
                                     af.sample_from_dirichlet(int(n_genes))]

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

                    gene.length = int(af.obtain_value(self.parameters["GENE_LENGTH"]))

            if intergenic_sequences == True:
                chromosome.obtain_locations()

            genome.chromosomes.append(chromosome)

            if interactome == True:
                genome.create_interactome()

        self.initial_genome = copy.deepcopy(genome)
        self.initial_gene_families =  copy.deepcopy(self.all_gene_families)


        return genome

    def read_genome(self, genome_file:str, intergenic_sequences = False,
                    family_rates = False, interactome = False):
        """
        Create a genome with genes and intergenic regions specified by the given
        `genome_file` (.gff).

        Parameters
        ----------
        genome_file : str
            the filename (.gff) that annotates the genes in the genome. A
            single chromosome is assumed.
        intergenic_sequences : bool, optional
            [description], by default False
        family_rates : bool, optional
            [description], by default False
        interactome : bool, optional
            [description], by default False

        Returns
        -------
        Genome
            the newly constructed genome with gene and intergenic sizes having
            the lengths specified in `genome_file`.

        Notes
        -----
            Only makes a single circular chromosome at the moment.
        """
        genome = Genome()
        genome.species = "Root"
        time = 0

        chrom_len, gene_features = af.parse_GFF(genome_file)

            #Create a chromosome of the appropriate shape:
        shape = "C"
        chromosome: Union[LinearChromosome, CircularChromosome]
        if shape == "L":
            chromosome = LinearChromosome(chrom_len)
            chromosome.shape = "L"
        elif shape == "C":
            chromosome = CircularChromosome(num_nucleotides=chrom_len)
            chromosome.shape = "C"

            #Create the genes:
        prev_feature = None     #previous feature used to make a gene
        for feature in gene_features:
            if prev_feature and prev_feature.location.end > feature.location.start:
                print(f'WARNING: skipping the creation of overlapping gene '
                      f'{feature.id},\n\tas it overlaps with {prev_feature.id}')
            else:
                gene, gene_family = self.make_gene(feature, genome.species,
                                                   time, family_rates,
                                                   self.parameters["RATE_FILE"] != "False")
                chromosome.genes.append(gene)
                self.all_gene_families[str(self.gene_families_counter)] = gene_family

                prev_feature = feature

            #Create the intergenes:
        if intergenic_sequences:    #first intergene is after the first gene
            chromosome.has_intergenes = True

            if shape == "L":        #before the first gene
                chromosome.intergenes.append(Intergene(chromosome.genes[0].start))

            for gene1, gene2 in af.pairwise(chromosome.genes):
                assert gene2.start >= gene1.end
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
            chromosome.obtain_locations()

        genome.chromosomes.append(chromosome)

        if interactome:
            genome.create_interactome()

        self.initial_genome = copy.deepcopy(genome)
        self.initial_gene_families =  copy.deepcopy(self.all_gene_families)

        return genome

    def make_gene(self, gene_feature: SeqFeature, species_tree_node: str,
                  time: int, family_mode = False, empirical_rates = False) -> Tuple[Gene, GeneFamily]:
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

        gene.gene_family = str(self.gene_families_counter)
        gene.species = species_tree_node
        gene.start = gene_feature.location.start
        gene.end = gene_feature.location.end
        gene.length = gene.end - gene.start

        gene_family = GeneFamily(gene_family_id, time)
        gene_family.length = gene.length
        gene_family.genes.append(gene)
        gene_family.gff_id = gene_feature.id
        gene.gene_id = gene_family.obtain_new_gene_id()

        self.all_gene_families[gene_family_id] = gene_family
        self.all_gene_families[gene.gene_family].register_event(str(time), "O", species_tree_node)

        if family_mode and not empirical_rates:
            d, t,l, _, _, _ = self.generate_new_rates()
            gene_family.rates["DUPLICATION"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        elif family_mode and empirical_rates:
            d, t, l = self.generate_empirical_rates()
            gene_family.rates["DUPLICATION"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        initial_gene = copy.deepcopy(gene)
        initial_gene.species = "Initial"
        gene_family.genes.append(initial_gene)

        return gene, gene_family

    def run(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSPOSITION"])

        o = af.obtain_value(self.parameters["ORIGINATION"])

        # First we prepare the first genome

        genome = self.fill_genome()

        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the original genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

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

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.evolve_genomes(d, t, l, i, c, o, current_time)

    def run_i(self):

        # Interactome mode

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSPOSITION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])
        rm = af.obtain_value(self.parameters["REMOVE"])
        rw = af.obtain_value(self.parameters["REWIRE"])

        # First we prepare the first genome

        genome = self.fill_genome(interactome = True)

        # We prepare to important dicts in this mode

        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the initial genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

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

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2


                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.evolve_genomes_i(d, t, l, i, c, o, rm, rw, current_time)

    def run_m(self):

        # First we prepare the first genome

        genome = self.fill_genome(family_rates=True)

        # We prepare to important dicts in this mode

        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the initial genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

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

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2


                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.evolve_genomes_m( current_time)


    def run_u(self):

        genome = self.fill_genome()
        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the original genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

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
                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.advanced_evolve_genomes(current_time)

    def run_f(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSPOSITION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])

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
        self.all_genomes["Root"] = genome

        # We add the original genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

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

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.advanced_evolve_genomes_f(d, t, l, i, c, o, current_time)

    def run_f_debug(self, injected_events): # Only for debugging purposes

        # First we prepare the root genome

        if self.root_genome_file:
            genome = self.read_genome(self.root_genome_file,
                                      intergenic_sequences=True)
        else:
            genome = self.fill_genome(intergenic_sequences=True)


        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the original genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

        all_species_tree_events = [("T",float(x[0]),x[1],x[2]) for x in self.tree_events]
        
        # Second, we compute the time to the next event:

        all_events = sorted(all_species_tree_events + injected_events, key=lambda x: x[1]) # We put together the two types of events and order by time

        for item in all_events:

            if item[0] == "T":
                etype, time, event, nodes = item

                if event == "S":

                    sp, c1, c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes

                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, time)

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(time)
                    break

            else:
                # If the event is a Genome level event
                etype, time, event, nodes, r = item
                c1, c2, d = r
                 
                if event == "D":

                    self.make_duplication_within_intergene(ch, c1, c2, d, lineage, time)

                elif event == "T":
                    pass

                elif event == "L":

                    pseudo = False
                    if numpy.random.uniform(0,1) <= float(self.parameters["PSEUDOGENIZATION"]):
                        pseudo = True
                    self.make_loss_intergenic(ch, c1, c2, d, lineage, time, pseudo)

                    return "L", lineage

                elif event == "I":

                    lineage = nodes
                    for ch in self.all_genomes[lineage]:
                        ch.obtain_flankings()
                        ch.obtain_locations()

                    #for intergene in ch.iter_intergenes():
                    #    print(intergene)

                    self.make_inversion_intergenic(ch, c1, c2, d, lineage, time)

                elif event == "P":

                    c3 = ch.select_random_intergenic_coordinate_excluding(c1, c2, d)
                    self.make_transposition_intergenic(ch, c1, c2, d, c3, lineage, time)

                    return "P", lineage

                elif event == "O":
                        
                    ch = self.all_genomes[lineage].select_random_chromosome()
                    intergene_coordinate = ch.select_random_coordinate_in_intergenic_regions()
                    self.make_origination_intergenic(ch, intergene_coordinate,
                                                     lineage, time)

                    return "O", lineage
    
    def generate_new_rates(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        p = af.obtain_value(self.parameters["TRANSPOSITION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])

        return d,t,l,i,p,o

    def generate_empirical_rates(self):

        mlen = len(self.empirical_rates)
        d,t,l = self.empirical_rates[numpy.random.randint(mlen)]

        return d,t,l

    def read_rates(self, rates_folder):

        self.branch_event_rates = dict()
        self.branch_extension_rates = dict()
        self.transfer_rates = dict()

        with open(os.path.join(rates_folder, "Event_rates.tsv")) as f:
            f.readline()
            for line in f:
                sp, d, t, l, i, c, o,  = line.split("\t")
                self.branch_event_rates[sp] = tuple([float(x) for x in (d, t, l, i, c, o)])

        with open(os.path.join(rates_folder, "Extension_rates.tsv")) as f:
            f.readline()
            for line in f:
                sp, d, t, l, i, c,  = line.split("\t")
                self.branch_extension_rates[sp] = tuple([x for x in (d, t, l, i, c)])

        with open(os.path.join(rates_folder, "Transfer_rates.tsv")) as f:
            f.readline()
            for line in f:
                dn, rc, wt  = line.split("\t")
                if dn not in self.transfer_rates:
                    self.transfer_rates[dn] = dict()
                if rc not in self.transfer_rates[dn]:
                    self.transfer_rates[dn][rc] = 0.0

                self.transfer_rates[dn][rc] = float(wt)


    def choose_event(self, duplication, transfer, loss, inversion, transposition, origination):

        draw = numpy.random.choice(["D", "T", "L", "I", "P", "O"], 1,
                                   p=af.normalize([duplication, transfer, loss, inversion, transposition, origination]))
        return draw

    def choose_event_i(self, duplication, transfer, loss, inversion, transposition, origination, remove, rewire):

        draw = numpy.random.choice(["D", "T", "L", "I", "P", "O", "RM", "RW"], 1,
                                   p=af.normalize([duplication, transfer, loss, inversion, transposition, origination, remove, rewire]))
        return draw

    def choose_recipient(self, lineages_alive, donor):
        possible_recipients = [x for x in lineages_alive if x != donor]
        if len(possible_recipients) > 1:
            recipient = random.choice(possible_recipients)
            return recipient
        else:
            return None

    def evolve_genomes(self, duplication, transfer, loss, inversion, transposition, origination, time):

        lineage = random.choice(list(self.active_genomes))
        event = self.choose_event(duplication, transfer, loss, inversion, transposition, origination)

        if event == "D":
            d_e = self.parameters["DUPLICATION_EXTENSION"]
            self.make_duplication(d_e, lineage, time)
            return "D", lineage

        elif event == "T":

            t_e = self.parameters["TRANSFER_EXTENSION"]

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                donor = lineage

                # We choose a recipient

                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = random.choice(possible_recipients)

                self.make_transfer(t_e, donor, recipient, time)
                return "T", donor+"->"+recipient

            else:
                return None

        elif event == "L":

            l_e = self.parameters["LOSS_EXTENSION"]

            self.make_loss(l_e, lineage, time)
            return "L", lineage

        elif event == "I":
            i_e = self.parameters["INVERSION_EXTENSION"]
            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "P":
            c_e = self.parameters["TRANSPOSITION_EXTENSION"]
            self.make_transposition(c_e, lineage, time)
            return "P",lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)
            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)
            return "O", lineage

    def evolve_genomes_i(self, duplication, transfer, loss, inversion, transposition, origination, remove, rewire, time):

        d_e = self.parameters["DUPLICATION_EXTENSION"]
        t_e = self.parameters["TRANSFER_EXTENSION"]
        l_e = self.parameters["LOSS_EXTENSION"]
        i_e = self.parameters["INVERSION_EXTENSION"]
        c_e = self.parameters["TRANSPOSITION_EXTENSION"]

        lineage = random.choice(list(self.active_genomes))
        event = self.choose_event_i(duplication, transfer, loss, inversion, transposition, origination, remove, rewire)

        if event == "D":

            self.make_duplication_interactome(d_e, lineage, time)
            return "D", lineage

        elif event == "T":

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                donor = lineage
                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = random.choice(possible_recipients)

                self.make_transfer_interactome(t_e, donor, recipient, time)
                return "T", donor + "->" + recipient

            else:
                return None

        elif event == "L":

            self.make_loss_interactome(l_e, lineage, time)
            return "L", lineage

        elif event == "I":
            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "P":
            self.make_transposition(c_e, lineage, time)
            return "P", lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)
            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)

            # We need to insert the gene in the interactome too, with preferential attachment

            interactome = self.all_genomes[lineage].interactome

            node_degrees = [d + 1 for n, d in interactome.degree()]
            choice = numpy.random.choice(interactome.nodes, 1, p=af.normalize(node_degrees))[0]

            interactome.add_node(str(gene))
            interactome.add_edge(str(gene), choice)

            return "O", lineage

        elif event == "RM":

            self.make_remove_edge(lineage, time)
            return "RM", lineage


        elif event == "RW":

            self.make_rewiring_edge(lineage, time)

            return "RW", lineage


    def evolve_genomes_m(self, time):

        d_e = self.parameters["DUPLICATION_EXTENSION"]
        t_e = self.parameters["TRANSFER_EXTENSION"]
        l_e = self.parameters["LOSS_EXTENSION"]
        i_e = self.parameters["INVERSION_EXTENSION"]
        c_e = self.parameters["TRANSPOSITION_EXTENSION"]

        ####


        ####

        mactive_genomes = list(self.active_genomes)
        mweights = list()
        for genome in mactive_genomes:
            lineage_weight = 0
            for chromosome in self.all_genomes[genome]:
                for gene in chromosome:
                    for r,vl in self.all_gene_families[gene.gene_family].rates.items():
                        lineage_weight += vl
            mweights.append(lineage_weight)

        lineage = numpy.random.choice(mactive_genomes, 1, p=af.normalize(mweights))[0]

        d, t, l, i, p, o = 0, 0, 0, 0, 0, 0

        for chromosome in self.all_genomes[lineage]:
            for gene in chromosome:
                d += self.all_gene_families[gene.gene_family].rates["DUPLICATION"]
                t += self.all_gene_families[gene.gene_family].rates["TRANSFER"]
                l += self.all_gene_families[gene.gene_family].rates["LOSS"]

        i += af.obtain_value((self.parameters["INVERSION"]))
        p += af.obtain_value((self.parameters["TRANSPOSITION"]))
        o += af.obtain_value((self.parameters["ORIGINATION"]))

        #print(d,t,l,i,p,o)

        event = self.choose_event(d, t, l, i, p, o)

        ####

        if event == "D":

            self.make_duplication(d_e, lineage, time, family_mode=True)
            return "D", lineage

        elif event == "T":

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                donor = lineage
                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = random.choice(possible_recipients)

                self.make_transfer(t_e, donor, recipient, time, family_mode = True)
                return "T", donor + "->" + recipient

            else:
                return None

        elif event == "L":
            self.make_loss(l_e, lineage, time, family_mode = True)
            return "L", lineage

        elif event == "I":
            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "P":
            self.make_transposition(c_e, lineage, time)
            return "P", lineage

        elif event == "O":

            if  self.parameters["RATE_FILE"] == "False":
                gene, gene_family = self.make_origination(lineage, time, family_mode=True)

            elif  self.parameters["RATE_FILE"] != "False":
                gene, gene_family = self.make_origination(lineage, time, family_mode=True,
                                                          empirical_rates=True)
            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)

            return "O", lineage


    def advanced_evolve_genomes(self, time):

        active_genomes = list(self.active_genomes)
        lineage = numpy.random.choice(active_genomes, 1, p=af.normalize(
            [sum(self.branch_event_rates[x]) for x in active_genomes]))[0]

        d,t,l,i,c,o = self.branch_event_rates[lineage]

        event = self.choose_event(d,t,l,i,c,o)

        d_e, t_e, l_e, i_e, c_e = self.branch_extension_rates[lineage]

        if event == "D":
            self.make_duplication(d_e, lineage, time)
            return "D", lineage

        elif event == "T":

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                recipient = self.choose_advanced_recipient(possible_recipients, lineage)
                if recipient != None:
                    donor = lineage
                    self.make_transfer(t_e, donor, recipient, time)
                    return "T", donor+"->"+recipient

            else:
                return None

        elif event == "L":

            self.make_loss(l_e, lineage, time)
            return "L", lineage

        elif event == "I":
            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "P":
            self.make_transposition(c_e, lineage, time)
            return "P",lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)

            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)

            return "O", lineage

    def advanced_evolve_genomes_f(self, duplication, transfer, loss, inversion, transposition, origination, time):
        
        # Evolve genomes with intergenes

        d_e = int(af.obtain_value(self.parameters["DUPLICATION_EXTENSION"]))
        t_e = int(af.obtain_value(self.parameters["TRANSFER_EXTENSION"]))
        l_e = int(af.obtain_value(self.parameters["LOSS_EXTENSION"]))
        i_e = int(af.obtain_value(self.parameters["INVERSION_EXTENSION"]))
        c_e = int(af.obtain_value(self.parameters["TRANSPOSITION_EXTENSION"]))

        distribution  = self.parameters["GENE_LENGTH"].split(":")[0]

        if distribution == "f" or distribution == "n":
            mean_gene_length = int(self.parameters["GENE_LENGTH"].split(":")[1].split(";")[0])
        elif distribution == "u":
            u1,u0 = self.parameters["GENE_LENGTH"].split(":")[1].split(";")
            mean_gene_length = (int(u1) - int(u0))/2
        else:
            print("Error, please switch the distribution type por the gene length")
            return 0

        mean_intergene_length = int(self.parameters["INTERGENE_LENGTH"])
        multiplier = 1.0 / mean_intergene_length


        lineage = random.choice(list(self.active_genomes))
        event = self.choose_event(duplication, transfer, loss, inversion, transposition, origination)

        for ch in self.all_genomes[lineage]:
            ch.obtain_flankings()
            ch.obtain_locations()

        if event == "D":

            r = self.select_advanced_length(lineage, 1/d_e * multiplier)
            if r == None:
                return None
            else:
                
                ch, c1, c2, d = r
                self.make_duplication_within_intergene(ch, c1, c2, d, lineage, time)

            return "D", lineage

        elif event == "T":


            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                donor = lineage
                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = random.choice(possible_recipients)

                r = self.select_advanced_length(lineage, 1/t_e * multiplier)

                if r == None:
                    return None

                ch, c1, c2, d = r

                chreceptor: CircularChromosome = self.all_genomes[recipient].select_random_chromosome()
                chreceptor.obtain_locations()
                c3 = chreceptor.select_random_coordinate_in_intergenic_regions()

                self.make_transfer_intergenic(ch, c1, c2, d, donor, chreceptor,
                                              c3, recipient, time)

                return "T", donor + "->" + recipient

            else:
                return None

        elif event == "L":

            r = self.select_advanced_length(lineage, 1/l_e * multiplier)

            if r == None:
                return None
            else:
                ch, c1, c2, d = r
                pseudo = False
                if numpy.random.uniform(0,1) <= float(self.parameters["PSEUDOGENIZATION"]):
                    pseudo = True
                self.make_loss_intergenic(ch, c1, c2, d, lineage, time, pseudo)

            return "L", lineage

        elif event == "I":

            r = self.select_advanced_length(lineage, 1/i_e * multiplier)

            if r == None:
                return None
            else:
                ch, c1, c2, d = r
                self.make_inversion_intergenic(ch, c1, c2, d, lineage, time)

            return "I", lineage

        elif event == "P":

            r = self.select_advanced_length(lineage, 1/c_e * multiplier)
            if r == None:
                return None

            ch, c1, c2, d = r

            c3 = ch.select_random_intergenic_coordinate_excluding(c1, c2, d)
            self.make_transposition_intergenic(ch, c1, c2, d, c3, lineage, time)

            return "P", lineage

        elif event == "O":
                
            ch = self.all_genomes[lineage].select_random_chromosome()
            intergene_coordinate = ch.select_random_coordinate_in_intergenic_regions()
            self.make_origination_intergenic(ch, intergene_coordinate,
                                             lineage, time)

            return "O", lineage


#    def advanced_evolve_genomes_f_debug(self, duplication, transfer, loss, inversion, transposition, origination, time):
#        
#        
#
#        for ch in self.all_genomes[lineage]:
#            ch.obtain_flankings()
#            ch.obtain_locations()
#
#        if event == "D":
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
#            return "D", lineage
#
#        elif event == "T":
#
#
#            # We choose a recipient
#
#            possible_recipients = [x for x in self.active_genomes if x != lineage]
#
#            if len(possible_recipients) > 0:
#
#                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
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
#                return "T", donor + "->" + recipient
#
#            else:
#                return None
#
#        elif event == "L":
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
#            return "L", lineage
#
#        elif event == "I":
#
#            r = self.select_advanced_length(lineage, 1/i_e * multiplier)
#
#            if r == None:
#                return None
#            else:
#                c1, c2, d = r
#                self.make_inversion_intergenic(ch, c1, c2, d, lineage, time)
#
#            return "I", lineage
#
#        elif event == "P":
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
#            return "P", lineage
#
#        elif event == "O":
#                
#            ch = self.all_genomes[lineage].select_random_chromosome()
#            intergene_coordinate = ch.select_random_coordinate_in_intergenic_regions()
#            self.make_origination_intergenic(ch, intergene_coordinate,
#                                             lineage, time)
#
#            return "O", lineage


    def get_time_to_next_event(self, n, events):

        total = 0.0
        for __ in range(n):
            total += sum(events)

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        else:
            time = numpy.random.exponential(1/total)
            return time

    def get_time_to_next_event_advanced_modes(self):
        # To obtain the time to next event in case that we have different rates per branch
        total = 0.0

        for lineage in self.active_genomes:
            total += sum(self.branch_event_rates[lineage])

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        time = numpy.random.exponential(1 / total)
        return time

    def get_time_to_next_event_family_mode(self):

        total = 0.0

        for lineage in self.active_genomes:
            for chromosome in self.all_genomes[lineage]:
                for gene in chromosome:
                    for r,vl in self.all_gene_families[gene.gene_family].rates.items():
                        total += vl

        total_active =  len(self.active_genomes)

        total += af.obtain_value((self.parameters["INVERSION"])) * total_active
        total += af.obtain_value((self.parameters["TRANSPOSITION"])) * total_active
        total += af.obtain_value((self.parameters["ORIGINATION"])) * total_active

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        time = numpy.random.exponential(1 / total)

        return time

    def increase_distances(self, time_to_next_event, active_lineages):

        for node in active_lineages:
            node.dist += time_to_next_event

    def make_origination(self, species_tree_node, time, family_mode = False, empirical_rates = False):

        self.gene_families_counter += 1
        gene_family_id = str(self.gene_families_counter)

        gene = Gene()
        gene.determine_orientation()

        gene.gene_family = str(self.gene_families_counter)
        gene.species = species_tree_node

        gene_family = GeneFamily(gene_family_id, time)
        gene_family.length = int(af.obtain_value(self.parameters["GENE_LENGTH"]))
        gene.length = gene_family.length

        gene_family.genes.append(gene)
        gene.gene_id = gene_family.obtain_new_gene_id()

        self.all_gene_families[gene_family_id] = gene_family
        self.all_gene_families[gene.gene_family].register_event(str(time), "O", species_tree_node)

        if family_mode == True and empirical_rates == False:

            d, t,l, _, _, _ = self.generate_new_rates()
            gene_family.rates["DUPLICATION"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        elif family_mode == True and empirical_rates == True:
            d, t, l = self.generate_empirical_rates()
            gene_family.rates["DUPLICATION"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        return gene, gene_family

    def make_origination_intergenic(self, chromosome: CircularChromosome, c: int,
                                    lineage, time) -> Gene:

        gene, _ = self.make_origination(lineage, time)
        location = chromosome.return_location_by_coordinate(c, within_intergene=True)
        chromosome.insert_gene_within_intergene(c, location, gene)

        orig = Origination(location, c, gene.length, lineage, time)
        chromosome.event_history.append(orig)

        return gene

    def make_speciation(self, sp, c1, c2, time, intergene=False):

        # This function receives a genome and the names of the two branching lineages of the species node

        genome_sp = self.all_genomes[sp]

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
                ch1 = CircularChromosome()
                ch2 = CircularChromosome()
                ch1.shape = "C"
                ch2.shape = "C"

            elif shape == "L":
                ch1 = LinearChromosome()
                ch2 = LinearChromosome()
                ch1.shape = "L"
                ch2.shape = "L"

            genome1.chromosomes.append(ch1)
            genome2.chromosomes.append(ch2)

            if chromosome.has_intergenes:
                ch1.has_intergenes = True
                ch2.has_intergenes = True

            for gene in chromosome:

                new_id1 = self.return_new_identifiers_for_segment([gene])
                new_id2 = self.return_new_identifiers_for_segment([gene])

                new_gene1 = af.copy_segment([Gene()], new_id1)[0]
                new_gene2 = af.copy_segment([Gene()], new_id2)[0]

                new_gene1.species = c1
                new_gene2.species = c2

                new_gene1.orientation = gene.orientation
                new_gene2.orientation = gene.orientation

                new_gene1.length = gene.length
                new_gene2.length = gene.length

                gene_family = self.all_gene_families[gene.gene_family]
                gene_family.genes.append(new_gene1)
                gene_family.genes.append(new_gene2)

                new_gene1.gene_family = gene.gene_family
                new_gene2.gene_family = gene.gene_family

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

                self.all_gene_families[gene.gene_family].register_event(str(time), "S", ";".join(map(str,nodes)))

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

    def make_extinction(self, sp, time, intergene=False):

        # We have to inactivate all the genes

        genome = self.all_genomes[sp]

        for chromosome in genome:
            for gene in chromosome:
                gene.active = False
                self.all_gene_families[gene.gene_family].register_event(str(time), "E", ";".join(map(str,[sp, gene.gene_id])))

        if intergene:
            pass

    def make_end(self, time):

        for genome_name in self.active_genomes:
            genome = self.all_genomes[genome_name]
            for chromosome in genome:
                for gene in chromosome:
                    gene.active = False
                    self.all_gene_families[gene.gene_family].register_event(str(time), "F", ";".join(
                        map(str, [genome.species, gene.gene_id])))


    def make_duplication(self, p, lineage, time, family_mode = False):

        chromosome = self.all_genomes[lineage].select_random_chromosome()

        if family_mode == True:
            affected_genes = chromosome.obtain_affected_genes_accounting_for_family_rates(p, self.all_gene_families, "DUPLICATION")
        else:
            affected_genes = chromosome.obtain_affected_genes(p)

        segment = chromosome.obtain_segment(affected_genes)

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        # We insert the two new segments after the last position of the old segment

        chromosome.insert_segment(affected_genes[-1], copied_segment1 + copied_segment2)

        # And we remove the old segment

        chromosome.remove_segment(segment)

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
            gene_family = gene.gene_family

            self.all_gene_families[gene_family].genes.append(copied_segment1[i])
            self.all_gene_families[gene_family].genes.append(copied_segment2[i])
            self.all_gene_families[gene_family].register_event(time, "D", ";".join(map(str, nodes)))

    def make_duplication_interactome(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        # We insert the two new segments after the last position of the old segment

        chromosome.insert_segment(affected_genes[-1], copied_segment1 + copied_segment2)

        # And we remove the old segment

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been a duplication

        # We need to create two dicts to change the names of the interactome

        new_genes_1 = dict()

        for i, gene in enumerate(segment):

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            gene.active = False

            gene_family = gene.gene_family

            self.all_gene_families[gene_family].genes.append(copied_segment1[i])
            self.all_gene_families[gene_family].genes.append(copied_segment2[i])

            self.all_gene_families[gene.gene_family].register_event(time, "D", ";".join(map(str, nodes)))

            # We update the new gene name (which by default is going to be the copied_segment1)

            gene1 = copied_segment1[i]
            gene2 = copied_segment2[i]

            new_genes_1[str(gene)] = str(gene1)

            self.all_genomes[lineage].interactome = nx.relabel_nodes(self.all_genomes[lineage].interactome,
                                                                     new_genes_1)

            # WE ADD THE NEW NODE

            self.all_genomes[lineage].interactome.add_node(str(gene2))

            # We distribute node depending on the parameter PROPORTION

            # If p is 1, all the links go to the first node
            # If p is 0, all the links go to the second node
            # If p is 0.5, equal repartition.

            PROPORTION = 0.5

            n_edges_to_old_node = int(PROPORTION * len(self.all_genomes[lineage].interactome.edges(str(gene1))))
            myedges = list(self.all_genomes[lineage].interactome.edges(str(gene1)))
            random.shuffle(myedges)
            edges_to_new_node = myedges[n_edges_to_old_node:]
            edges_to_add_to_new_node = [(str(gene2), x[1]) for x in edges_to_new_node]
            edges_to_remove_to_old_node = edges_to_new_node

            # Now, I have to remove the edges

            self.all_genomes[lineage].interactome.remove_edges_from(edges_to_remove_to_old_node)
            self.all_genomes[lineage].interactome.add_edges_from(edges_to_add_to_new_node)


    def make_duplication_within_intergene(self, chromosome: CircularChromosome,
                                          c1: int, c2: int, d: T_DIR,
                                          lineage: str, time: float):
        """
        Do a duplication that acts on the given pair of intergene spcific
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
        
        chromosome: CircularChromosome = self.all_genomes[lineage].select_random_chromosome()
        r = chromosome.return_affected_region(c1, c2, d)

        if r == None:
            return None

        genepositions, intergenepositions, leftlengths, rightlengths, int1, int2 = r
        segment = chromosome.obtain_segment(genepositions)
        intergene_segment = chromosome.obtain_intergenic_segment(intergenepositions[1:])

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # We duplicate the genes

        new_segment_1 = af.copy_segment(segment, new_identifiers1)
        new_segment_2 = af.copy_segment(segment, new_identifiers2)

        # And the intergenes, which will have incorrect lengths until we update
        # them later

        new_intergene_segment_1 = [copy.deepcopy(chromosome.intergenes[x]) for x in intergenepositions[1:]]
        new_intergene_segment_2 = [copy.deepcopy(chromosome.intergenes[x]) for x in intergenepositions[1:]]

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

        if d == LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        dup = TandemDup(int1, int2, c1, c2, len(new_intergene_segment_1),
                        specificlen, totallen, lineage, time)
            
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

            gene_family = gene.gene_family

            self.all_gene_families[gene_family].genes.append(new_segment_1[i])
            self.all_gene_families[gene_family].genes.append(new_segment_2[i])

            self.all_gene_families[gene.gene_family].register_event(time, "D", ";".join(map(str, nodes)))

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
        chromosome.obtain_flankings()

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
            ca = self.complete_tree.get_common_ancestor(mrecipient, mdonor).name
            x1 = self.distances_to_start[ca]
            td =  time - x1
            weights.append(td)

        beta = min(alpha * af.normalize(weights))
        val = (alpha * af.normalize(weights)) - beta
        pvector = af.normalize(numpy.exp(-val))

        draw = numpy.random.choice(possible_recipients, 1, p=pvector)[0]

        return draw

    def choose_advanced_recipient(self, possible_recipients, donor):

        weights = list()

        for recipient in possible_recipients:
            weights.append(self.transfer_rates[donor][recipient])

        if sum(weights) == 0:
            return None


        draw = numpy.random.choice(possible_recipients, 1, p=af.normalize(weights))[0]

        return draw

    def make_transfer(self, p, donor, recipient, time, family_mode = False):

        chromosome1 = self.all_genomes[donor].select_random_chromosome()

        if family_mode == True:
            affected_genes = chromosome1.obtain_affected_genes_accounting_for_family_rates(p, self.all_gene_families, "TRANSFER")
        else:
            affected_genes = chromosome1.obtain_affected_genes(p)

        segment = chromosome1.obtain_segment(affected_genes)
        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        inverted = False

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        # We insert the first segment (leaving transfer) in the same position than the previous segment
        # We do this just to change the identifiers of the numbers

        chromosome1.insert_segment(affected_genes[0], copied_segment1)

        # And we remove the old segment

        chromosome1.remove_segment(segment)

        # Now we insert the transfer segment in the recipient genome in one of the homologous position.

        if numpy.random.uniform(0,1) <= self.parameters["REPLACEMENT_TRANSFER"]:

            possible_positions = list()

            for chromosome in self.all_genomes[recipient]:
                for direction, positions in chromosome.get_homologous_position(segment):
                    possible_positions.append((direction, positions, chromosome))

            if len(possible_positions) != 0:

                direction, positions, chromosome2 = random.choice(possible_positions)


                if direction == "F":

                    # I replace gene by gene the segment

                    i = 0

                    for position in positions:

                        # First I inactivate the gene

                        gene = chromosome2.genes[position]
                        gene.active = False
                        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(
                            map(str, [recipient, gene.gene_id])))

                        # And then I replace

                        chromosome2.genes[position] = copied_segment2[i]

                        i += 1

                elif direction == "B":
                    # I invert the segment and I replace gene by gene the segment
                    inverted = True

                    copied_segment2 = copied_segment2[::-1]

                    for gene in copied_segment2:
                        gene.change_sense()

                    i = 0
                    for position in positions:
                        # First I inactivate the gene
                        gene = chromosome2.genes[position]
                        gene.active = False
                        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(
                            map(str, [recipient, gene.gene_id])))

                        # And then I replace

                        chromosome2.genes[position] = copied_segment2[i]


                        i += 1
            else:

                # Normal transfers
                chromosome2 = self.all_genomes[recipient].select_random_chromosome()
                position = chromosome2.select_random_position()
                chromosome2.insert_segment(position, copied_segment2)
        else:
            # Normal transfer
            chromosome2 = self.all_genomes[recipient].select_random_chromosome()
            position = chromosome2.select_random_position()
            chromosome2.insert_segment(position, copied_segment2)

        # We have to register in the affected gene families that there has been a transfer event

        if inverted == True:
            # We invert again to store the event
            copied_segment2 = copied_segment2[::-1]

        for i, gene in enumerate(segment):

            gene.active = False

            # The code for the node is:
            # 1. Branch of the species tree for the donor genome
            # 2. Id of the gene that is transferred
            # 3. Id of the gene that remains in the donor genome
            # 4. Branch of the species tree for the recipient genome
            # 5. Id of the new gene arriving

            copied_segment1[i].species = donor
            copied_segment2[i].species = recipient

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            self.all_gene_families[gene.gene_family].register_event(time, "T", ";".join(map(str,nodes)))

    def make_transfer_interactome(self, p, donor, recipient, time):

        chromosome1 = self.all_genomes[donor].select_random_chromosome()
        affected_genes = chromosome1.obtain_affected_genes(p)
        segment = chromosome1.obtain_segment(affected_genes)

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        inverted = False

        is_replacement_transfer = False
        replaced_genes = list()

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        # We insert the first segment (leaving transfer) in the same position than the previous segment
        # We do this just to change the identifiers of the numbers

        chromosome1.insert_segment(affected_genes[0], copied_segment1)

        # And we remove the old segment

        chromosome1.remove_segment(segment)

        # Now we insert the transfer segment in the recipient genome in one of the homologous position.

        if numpy.random.uniform(0,1) <= self.parameters["REPLACEMENT_TRANSFER"]:

            possible_positions = list()

            for chromosome in self.all_genomes[recipient]:
                for direction, positions in chromosome.get_homologous_position(segment):
                    possible_positions.append((direction, positions, chromosome))

            if len(possible_positions) != 0:

                direction, positions, chromosome2 = random.choice(possible_positions)


                if direction == "F":

                    # I replace gene by gene the segment

                    i = 0

                    for position in positions:

                        # First I inactivate the gene

                        gene = chromosome2.genes[position]
                        gene.active = False


                        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(
                            map(str, [recipient, gene.gene_id])))
                        # And then I replace
                        chromosome2.genes[position] = copied_segment2[i]
                        i += 1
                        is_replacement_transfer = True

                        replaced_genes.append((gene, copied_segment2[i]))

                elif direction == "B":
                    # I invert the segment and I replace gene by gene the segment
                    inverted = True
                    copied_segment2 = copied_segment2[::-1]
                    for gene in copied_segment2:
                        gene.change_sense()

                    i = 0
                    for position in positions:
                        # First I inactivate the gene
                        gene = chromosome2.genes[position]
                        gene.active = False

                        replaced_genes.append(gene)

                        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(
                            map(str, [recipient, gene.gene_id])))
                        # And then I replace
                        chromosome2.genes[position] = copied_segment2[i]
                        i += 1
                        is_replacement_transfer = True

                        replaced_genes.append((gene, copied_segment2[i]))

            else:

                # Normal transfers
                chromosome2 = self.all_genomes[recipient].select_random_chromosome()
                position = chromosome2.select_random_position()
                chromosome2.insert_segment(position, copied_segment2)
        else:
            # Normal transfer
            chromosome2 = self.all_genomes[recipient].select_random_chromosome()
            position = chromosome2.select_random_position()
            chromosome2.insert_segment(position, copied_segment2)

        # We have to register in the affected gene families that there has been a transfer event

        if inverted == True:
            # We invert again to store the event
            copied_segment2 = copied_segment2[::-1]


        ## We register the event

        new_names_1 = dict()
        new_names_2 = dict()

        for i, gene in enumerate(segment):

            gene.active = False

            # The code for the node is:
            # 1. Branch of the species tree for the donor genome
            # 2. Id of the gene that is transferred
            # 3. Id of the gene that remains in the donor genome
            # 4. Branch of the species tree for the recipient genome
            # 5. Id of the new gene arriving

            copied_segment1[i].species = donor
            copied_segment2[i].species = recipient

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            self.all_gene_families[gene.gene_family].register_event(time, "T", ";".join(map(str,nodes)))

            new_names_1[str(gene)] = str(copied_segment1[i])

            ## We update the interactome

            # First we update the interactome in the donor lineage

            self.all_genomes[donor].interactome = nx.relabel_nodes(self.all_genomes[donor].interactome, new_names_1)

            # Second we update the interactome in the recipient lineage

            if is_replacement_transfer == True:

                ## All links pass now to the new gene
                ## The old gene gets no links

                new_names_2 = {str(n1):str(n2) for n1,n2 in replaced_genes}
                #new_names_2[str(gene)] = str(copied_segment2[i])
                self.all_genomes[recipient].interactome = nx.relabel_nodes(self.all_genomes[recipient].interactome, new_names_2)


            else:

                # It is not a replacement transfer. Preferential attachment

                node_degrees = [d + 1 for n, d in self.all_genomes[recipient].interactome.degree()]
                choice = numpy.random.choice(self.all_genomes[recipient].interactome.nodes, 1, p=af.normalize(node_degrees))[0]
                self.all_genomes[recipient].interactome.add_node(str(copied_segment2[i]))
                self.all_genomes[recipient].interactome.add_edge(str(copied_segment2[i]), choice)



    def make_transfer_intergenic(self, donorchrom: Chromosome,
                                 c1: int, c2: int, d: T_DIR, donor: str,
                                 receptorchrom: Chromosome, c3: int,
                                 receptor: str, time: int):

        r = donorchrom.return_affected_region(c1, c2, d)

        if r == None:
            return None

        else:
            gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r
            segment = donorchrom.obtain_segment(gpositions)

            # Get lengths from last intergene.
        specificlen = donorchrom.intergenes[-1].specific_flanking[1]
        totallen = donorchrom.intergenes[-1].total_flanking[1]
        numintergenes = len(donorchrom.intergenes)
            
        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        new_intergene_segment = [copy.deepcopy(donorchrom.intergenes[x])
                                 for x in igpositions[1:]]

        # We insert the first segment (leaving transfer) in the same position as the previous segment
        # We do this just to change the identifiers of the numbers

        # We insert in the same place

        position = gpositions[-1] + 1

        for i, gene in enumerate(copied_segment1):
            donorchrom.genes.insert(position + i, gene)

        # We remove the old copies:

        donorchrom.remove_segment(segment)

        # Normal transfer

        int3 = receptorchrom.return_location_by_coordinate(c3, True)
        position = int3.position + 1

        for i, gene in enumerate(copied_segment2):
            receptorchrom.genes.insert(position + i, gene)
        for i, intergene in enumerate(new_intergene_segment):
            receptorchrom.intergenes.insert(position + i, intergene)

        cut_position = (c3 - int3.sc1, int3.sc2 - c3)

        scar1 = receptorchrom.intergenes[int3.position]
        scar2 = receptorchrom.intergenes[position + i]

        if d == LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        scar1.length = leftlengths[1] + cut_position[0]
        scar2.length = rightlengths[0] + cut_position[1]

        tr = Transfer(int1, int2, c1, c2, specificlen, totallen, donor,
                      numintergenes, int3, c3, receptor, time)
        receptorchrom.event_history.append(tr)

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

            self.all_gene_families[gene.gene_family].register_event(time, "T", ";".join(map(str, nodes)))


    def make_loss(self, p, lineage, time, family_mode = False):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        if family_mode == True:
            affected_genes = chromosome.obtain_affected_genes_accounting_for_family_rates(p, self.all_gene_families, "LOSS")
        else:
            affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)

        # Now we check we are not under the minimum size

        if len(chromosome) - len(affected_genes) <= self.parameters["MIN_GENOME_SIZE"]:
            return 0

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned

        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))

    def make_loss_intergenic(self, chromosome, c1, c2, d: T_DIR, lineage, time,
                             pseudo=False):

        r = chromosome.return_affected_region(c1, c2, d)

        if r == None:
            return None

        gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r

        segment = chromosome.obtain_segment(gpositions)
        intergene_segment = chromosome.obtain_intergenic_segment(igpositions[1:])

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

        if d == LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        pseudo_intergenes = []
        if pseudo:
            pseudo_intergenes = intergene_segment

        loss = Loss(int1, int2, c1, c2, specificlen, totallen, lineage, time,
                    pseudo, pseudo_intergenes)
        chromosome.event_history.append(loss)

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
            self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))

    #def make_loss_family_mode(self, p, lineage, time):

    #    chromosome: CircularChromosome = self.all_genomes[lineage].select_random_chromosome()

    #    affected_genes = chromosome.obtain_affected_genes_accounting_for_family_rates(p, interactome)
    #    segment = chromosome.obtain_segment(affected_genes)

    #    # Now we check we are not under the minimum size

    #    if len(chromosome) - len(affected_genes) <= self.parameters["MIN_GENOME_SIZE"]:
    #        return 0

    #    chromosome.remove_segment(segment)

    #    # We have to register in the affected gene families that there has been as loss
    #    # All genes affected must be returned

    #    for gene in segment:
    #        gene.active = False
    #        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))

    def make_loss_interactome(self, p, lineage, time):

        interactome = self.all_genomes[lineage].interactome
        chromosome = self.all_genomes[lineage].select_random_chromosome()

        affected_genes = chromosome.obtain_affected_genes_accounting_for_connectedness(p, interactome)
        segment = chromosome.obtain_segment(affected_genes)

        # Now we check we are not under the minimum size

        if len(chromosome) - len(affected_genes) <= 0:
            return 0

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned

        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))
            # We remove from the connectome

            interactome.remove_node(str(gene))

    def make_inversion(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)
        chromosome.invert_segment(affected_genes)

        for gene in segment:
            self.all_gene_families[gene.gene_family].register_event(str(time), "I", ";".join(map(str,[lineage, gene.gene_id])))

    
    def make_inversion_intergenic(self, chromosome: CircularChromosome, c1: int,
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

        Parameters
        ----------
        c1 : int
            the first intergene specific breakpoint coordinate
        c2 : int
            the second intergene specific breakpoint coordinate
        d : T_DIR
            the direction, either left or right
        lineage : str
            the linege, which is name of the pendnt node
        time : float
            the time stamp of the event
        """
        
        r = chromosome.return_affected_region(c1, c2, d)

        if r == None:
            return None

        gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r

            # Get lengths from last intergene before modifying chromosome.
        specificlen = chromosome.intergenes[-1].specific_flanking[1]
        totallen = chromosome.intergenes[-1].total_flanking[1]

        segment = chromosome.obtain_segment(gpositions)
        chromosome.invert_segment(gpositions)

        if d == LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        sleftlen, srightlen, tleftlen, trightlen = 0, 0, 0, 0
        if gpositions[0] > gpositions[-1]:      #The inversion wraps:
            sleftlen, srightlen, tleftlen, trightlen = \
                 chromosome.inversion_wrap_lengths(gpositions)

        inv = Inversion(int1, int2, c1, c2, specificlen, totallen, sleftlen,
                        tleftlen, srightlen, trightlen, lineage, time)

        chromosome.event_history.append(inv)

        scar1 = chromosome.intergenes[igpositions[0]]
        scar2 = chromosome.intergenes[igpositions[-1]]

        scar1.length = leftlengths[0] + rightlengths[0]
        scar2.length = rightlengths[1] + leftlengths[1]

        assert scar1.length == len(inv.afterL)
        assert scar2.length == len(inv.afterR)      

        for gene in segment:
            self.all_gene_families[gene.gene_family].register_event(str(time), "I", ";".join(map(str,[lineage, gene.gene_id])))

    def make_transposition(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)
        chromosome.cut_and_paste(affected_genes)

        for gene in segment:
            self.all_gene_families[gene.gene_family].register_event(str(time), "P", ";".join(map(str,[lineage, gene.gene_id])))

    def make_transposition_intergenic(self, chromosome, c1, c2, d: T_DIR, c3,
                                      lineage, time):

        r = chromosome.return_affected_region(c1, c2, d)

        if r == None:
            return None

        gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r

        segment = chromosome.obtain_segment(gpositions)
        intergene_segment = chromosome.obtain_intergenic_segment(igpositions[1:])

        interval3 = chromosome.return_location_by_coordinate(c3, True)

        scar1 = chromosome.intergenes[igpositions[0]]
        scar2 = chromosome.intergenes[interval3.position]
        scar3 = chromosome.intergenes[igpositions[-1]]

        new_segment = list()
        new_intergene_segment = list()

        # Get old lengths from last intergene before modifying chromosome.

        specificlen = chromosome.intergenes[-1].specific_flanking[1]
        totallen = chromosome.intergenes[-1].total_flanking[1]
        numintergenes = len(chromosome.intergenes)

        # If we insert in the intergene i, the gene must occupy the position i - 1
        # We store it for reference

        left_gene = chromosome.genes[interval3.position]

        # Now we pop the genes

        for gene in segment:
            new_segment.append(chromosome.genes.pop(chromosome.genes.index(gene)))

        # And now we insert the genes at the right of the gene we saved before

        position = chromosome.genes.index(left_gene) + 1

        for i, gene in enumerate(new_segment):
            chromosome.genes.insert(position + i, gene)

        # We move the intergene on the right also

        # We save the position for insertion

        left_intergene = chromosome.intergenes[interval3.position]

        for intergene in intergene_segment:
            new_intergene_segment.append(chromosome.intergenes.pop(chromosome.intergenes.index(intergene)))

        # And now we insert the genes at the right of the gene we saved before

        position = chromosome.intergenes.index(left_intergene) + 1

        for i, intergene in enumerate(new_intergene_segment):
            chromosome.intergenes.insert(position + i, intergene)

        # Finally, we modify the segments so that they have the right length

        r5 = (c3 - interval3.sc1, interval3.sc2 - c3)

        if d == LEFT:
            leftlengths, rightlengths = rightlengths, leftlengths
            int1, int2 = int2, int1
            c1, c2 = c2, c1

        scar1.length = leftlengths[0] + rightlengths[1]
        scar2.length = leftlengths[1] + r5[0]
        scar3.length = rightlengths[0] + r5[1]

        trans = Transposition(int1, int2, c1, c2, interval3, c3, numintergenes,
                              specificlen, totallen, lineage, time)

        chromosome.event_history.append(trans)

        assert len(trans.afterH) == scar1.length
        assert len(trans.afterL) == scar2.length
        assert len(trans.afterR) == scar3.length

        for i, gene in enumerate(segment):
            self.all_gene_families[gene.gene_family].register_event(str(time), "P", ";".join(map(str,[lineage, gene.gene_id])))


    def make_rewiring_edge(self, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        position = chromosome.select_random_position()
        interactome = self.all_genomes[lineage].interactome
        normalized_weights = af.normalize([d + 1 for n, d in interactome.degree()])
        n1 = chromosome.genes[position]
        n2 = numpy.random.choice(interactome.nodes, 1, p=normalized_weights)[0]

        while (str(n1) == str(n2)):
            n2 = numpy.random.choice(interactome.nodes, 1, p=normalized_weights)[0]

        self.all_genomes[lineage].interactome.add_edge(str(n1), str(n2))
        self.all_gene_families[n1.gene_family].register_event(str(time), "RW", ";".join(map(str, [lineage, n1, n2])))

    def make_remove_edge(self, lineage, time):

        myedges = list(self.all_genomes[lineage].interactome.edges())

        if len(myedges) == 0:
            return None

        myedge = myedges[random.randint(0, len(myedges) - 1)]
        self.all_genomes[lineage].interactome.remove_edge(*myedge)

        self.all_gene_families[myedge[0].split("_")[0]].register_event(str(time), "RM", ";".join([lineage, myedge[0], myedge[1]]))


    def get_gene_family_tree(self):

        if len(self.gene_family["Gene_tree"].get_leaves()) < 3:
            return "None"
        else:
            return self.gene_family["Gene_tree"].write(format=1)


    def select_advanced_length(self, lineage: str, p: float, reps=100,
                               ) -> Tuple[CircularChromosome, int, int, T_DIR]:
        """
        Return a pair of specific coordinates for intergenic regions according
        to `p` on a chromosome of `lineage`.

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
        chromosome: CircularChromosome = self.all_genomes[lineage].select_random_chromosome()
            #The total number of intergenic nucleotides can be retrieved from
            #the last intergenic location:
        assert chromosome.map_of_locations[-1].isIntergenic()
        intergenic_specific_length = chromosome.map_of_locations[-1].sc2

        success = False
        counter = 0
        while counter <= reps and success == False:
            counter += 1

            sc1 = chromosome.select_random_coordinate_in_intergenic_regions()
            d = numpy.random.choice((LEFT, RIGHT), p=[0.5, 0.5])

            extension = numpy.random.geometric(p)

            if d == RIGHT:

                if sc1 + extension > intergenic_specific_length:
                    sc2 = sc1 + extension - intergenic_specific_length - 1
                    if sc2 < sc1:       # The event wraps to the right and
                        success = True  # doesn't cover the whole genome
                else:
                    sc2 = sc1 + extension
                    success = True

            elif d == LEFT:

                if sc1 - extension < 0:
                    sc2 = intergenic_specific_length - (extension - sc1) + 1
                    if sc1 < sc2:       # The event wraps to the left and
                        success = True  # doesn't cover the whole genome
                else:
                    sc2 = sc1 - extension
                    success = True

            if success:
                l1 = chromosome.return_location_by_coordinate(sc1, True)
                l2 = chromosome.return_location_by_coordinate(sc2, True)
                if l1 != l2:
                    return chromosome, sc1, sc2, d

                success = False

        return None

    def obtain_divisions(self):
        """
        Obtain the divisions at the root
        """
        # (Need to traverse in a different way once we simulate transfers)
        for node in self.complete_tree.traverse("postorder"):
            
            genome = self.all_genomes[node.name]            

            node.add_feature("cuts", set())
            
            for chromosome in genome:
                
                chromosome.obtain_flankings()

                if not node.is_leaf():
                    n1,n2 = node.get_children()
                    node.cuts = node.cuts.union(n1.cuts, n2.cuts)

                new_cuts = set()
                
                # This maps to the top of the branch all the cuts along that
                # branch.  We traverse from the end to the beginning:
                reverse_history = list(reversed(chromosome.event_history))
                for i, event1 in enumerate(reverse_history):
                    sc1 = event1.sbpL
                    sc2 = event1.sbpR
                    print("Mapping back in time these events")
                    print(sc1, sc2)

                    # Map the cuts through the events that are above:
                    for event2 in reverse_history[i+1:]:
                        sc1 = event2.afterToBeforeS(sc1)
                        sc2 = event2.afterToBeforeS(sc2)
                        print(sc1, sc2)
                        
                    new_cuts.add(sc1)
                    new_cuts.add(sc2)
                print("Done!")

                # This maps to the top of the branch all the cuts that were
                # acquired from previous nodes
                for cut in node.cuts:
                    for event in reverse_history:
                        cut = event.afterToBeforeS(cut)

                    new_cuts.add(cut)

                node.cuts = new_cuts
                  
        initial_chromosome = self.initial_genome.chromosomes[0]   
        all_cuts = list()
        self.natural_cuts = list()

        # These are the natural cuts from the intergenes (the limits with the
        # genes):
        for intergene in initial_chromosome.iter_intergenes():
            
            cut1, cut2 = intergene.specific_flanking
            self.natural_cuts.append((cut1, cut2))
            
            all_cuts.append(cut1)
            all_cuts.append(cut2)
        
        print("Natural cuts")
        for x in self.natural_cuts:
            print(x)
        print("***")
        print("Events cuts")
        print(sorted(list(self.complete_tree.cuts)))
        print("***")
        # These are the cuts from the events ## NOTE: NEED TO FIX THIS        
        all_cuts +=  sorted(list(self.complete_tree.cuts))
        
        # These are the cuts surrounding the Genes. We don't need to treat these
        cuts_to_ignore = {(x1[1], x2[0]) for x1, x2 in zip(self.natural_cuts,
                                                           self.natural_cuts[1:] +
                                                           [self.natural_cuts[0]])}
        
       
        all_cuts =  sorted(list(set(all_cuts))) # To remove possible repeated values
        initial_specific_flankings = zip(all_cuts, all_cuts[1:] + [all_cuts[0]])
        fam_id = 0
        
        self.initial_divisions = list()         # For debugging purposes

        for initial_specific_flanking in initial_specific_flankings:
            c1, c2 = initial_specific_flanking

            # We need to ignore the cuts where c1 is the right most extreme of
            # an intergene and c2 is the left most extreme of the next intergene
            if (c1, c2) in cuts_to_ignore:
                continue
            self.initial_divisions.append((c1,c2))
            fam_id += 1
            intergene = initial_chromosome.return_intergene_by_coordinate(c1)
            division = intergene.create_division("1", fam_id, initial_specific_flanking) # Identifier is 1 in the beginning
            division_family = DivisionFamily(fam_id, initial_specific_flanking)
            division_family.register_event(0, "O", "Root") # We register the origination 
            
            self.all_division_families[str(fam_id)] = division_family

    def obtain_events_for_divisions(self):

        """
        This function obtain all the events (at the Species level and at the Genome level) that occur in every division
        Once we have all the events for every division, we can reconstruct the division trees
        The simulation is run once more from zero, using the same events generated in the main simulation
        
        """

        def make_speciation_divisions(time, lineages):

            pn, c1, c2 = lineages.split(";") # Parent, child1, child2
            
            genome_pn = self.all_genomes_second[pn]
            genome1 = Genome()
            genome2 = Genome()
            
            self.all_genomes_second[c1] = genome1
            self.all_genomes_second[c2] = genome2

            for chromosome in genome_pn:

                ch1 = CircularChromosome()
                ch2 = CircularChromosome()

                genome1.chromosomes.append(ch1)
                genome2.chromosomes.append(ch2)

                ch1.has_intergenes = True
                ch2.has_intergenes = True

                for gene in chromosome:

                    new_id1 = self.return_new_identifiers_for_segment_with_divisions([gene]) 
                    new_id2 = self.return_new_identifiers_for_segment_with_divisions([gene])

                    new_gene1 = af.copy_segment([Gene()], new_id1)[0]
                    new_gene2 = af.copy_segment([Gene()], new_id2)[0]

                    new_gene1.species = c1
                    new_gene2.species = c2

                    new_gene1.orientation = gene.orientation
                    new_gene2.orientation = gene.orientation

                    new_gene1.length = gene.length
                    new_gene2.length = gene.length

                    gene_family = self.all_gene_families[gene.gene_family]
                    gene_family.genes.append(new_gene1)
                    gene_family.genes.append(new_gene2)

                    new_gene1.gene_family = gene.gene_family
                    new_gene2.gene_family = gene.gene_family

                    ch1.genes.append(new_gene1)
                    ch2.genes.append(new_gene2)

                for intergene in chromosome.iter_intergenes():
                    
                    intergene1 = copy.deepcopy(intergene)
                    intergene2 = copy.deepcopy(intergene)
                    
                    ch1.intergenes.append(intergene1)
                    ch2.intergenes.append(intergene2)

                    for division, division1, division2 in zip(intergene, intergene1, intergene2):

                        division1.identity = self.all_division_families[str(division1.division_family)].obtain_new_identifier()
                        division2.identity = self.all_division_families[str(division2.division_family)].obtain_new_identifier()

                        nodes = [pn, 
                                 division.identity,
                                 c1,
                                 division1.identity,
                                 c2,
                                 division2.identity
                                 ]
                        self.all_division_families[str(division.division_family)].register_event(str(time), "S", ";".join(map(str,nodes)))

            genome1.update_genome_species(c1)
            genome2.update_genome_species(c2)  

            #print(len([x for x in ch1.iter_intergenes()]))
            

        def make_extinction_divisions(time, lineage):
            genome = self.all_genomes_second[lineage]
            for ch in genome:
                for intergene in ch.iter_intergenes():
                    for division in intergene:
                        self.all_division_families[str(division.division_family)].register_event(str(time), "E", ";".join(
                        map(str, [genome.species, division.identity])))
        
        def make_end_divisions(time, lineage):
            
            genome = self.all_genomes_second[lineage]
            for ch in genome:
                for intergene in ch.iter_intergenes():
                    for division in intergene:
                        
                        self.all_division_families[str(division.division_family)].register_event(str(time), "F", ";".join(
                        map(str, [genome.species, division.identity])))

        def make_duplication_divisions(time, event):
            
            lineage = event.lineage
            chromosome = [x for x in self.all_genomes_second[lineage]][0] # This should be corrected, only works if the genomes have a single chromosome
            
            chromosome.obtain_flankings()
            chromosome.obtain_locations()

            c1 = event.sbpL
            c2 = event.sbpR

            r = chromosome.return_affected_region(event.sbpL, event.sbpR, RIGHT)

            if r == None:
                return None
    
            genepositions, intergenepositions, leftlengths, rightlengths, int1, int2 = r
            segment = chromosome.obtain_segment(genepositions)
            #intergene_segment = chromosome.obtain_intergenic_segment(intergenepositions[1:])

            new_identifiers1 = self.return_new_identifiers_for_segment_with_divisions(segment)
            new_identifiers2 = self.return_new_identifiers_for_segment_with_divisions(segment)

            # We duplicate the genes

            new_segment_1 = af.copy_segment(segment, new_identifiers1)
            new_segment_2 = af.copy_segment(segment, new_identifiers2)

            # And the intergenes, which will have incorrect lengths until we update
            # them later
            
            old_intergene_segment = [chromosome.intergenes[x] for x in intergenepositions[1:]]
            new_intergene_segment_1 = [copy.deepcopy(chromosome.intergenes[x]) for x in intergenepositions[1:]]
            new_intergene_segment_2 = [copy.deepcopy(chromosome.intergenes[x]) for x in intergenepositions[1:]]

            for intergene, intergene1, intergene2 in zip(old_intergene_segment, new_intergene_segment_1, new_intergene_segment_2):
                for division, division1, division2 in zip(intergene, intergene1, intergene2):
                    division1.identity = self.all_division_families[str(division1.division_family)].obtain_new_identifier()
                    division2.identity = self.all_division_families[str(division2.division_family)].obtain_new_identifier()
                    nodes = [lineage, 
                            division.identity,
                            lineage,
                            division1.identity,
                            lineage,
                            division2.identity
                            ]
                    self.all_division_families[str(division.division_family)].register_event(str(time), "D", ";".join(map(str,nodes)))


            scar0 = chromosome.intergenes[intergenepositions[0]]
            scar1 = new_intergene_segment_1[-1]
            scar2 = new_intergene_segment_2[-1]

            # Get old lengths from last intergene before modifying chromosome.
            
            specificlen = chromosome.intergenes[-1].specific_flanking[1]
            totallen = chromosome.intergenes[-1].total_flanking[1]

            #Replace the original set of genes/intergenes with the first copies:

            for pos, newgene in zip(genepositions, new_segment_1):
                chromosome.genes[pos] = newgene
            for pos, newintergene in zip(intergenepositions[1:], new_intergene_segment_1): # I need to obtain new identifiers for the divisions
                chromosome.intergenes[pos] = newintergene
                
            position = genepositions[-1] + 1

            # Insert the second copy of new genes and intergenes after the first.
            
            for i, gene in enumerate(new_segment_2):
                chromosome.genes.insert(position + i, gene)
            for i, intergene in enumerate(new_intergene_segment_2):
                chromosome.intergenes.insert(position + i, intergene)

            dup = TandemDup(int1, int2, c1, c2, len(new_intergene_segment_1),
                            specificlen, totallen, lineage, time)
                
            scar1.length = leftlengths[1] + rightlengths[0]
            scar2.length = rightlengths[1] + rightlengths[0]

            assert scar1.length == len(dup.afterC)
            assert scar2.length == len(dup.afterR)
            
            chromosome.event_history.append(dup)

            for i, gene in enumerate(segment):
                nodes = [gene.species,
                        gene.gene_id,
                        new_segment_1[i].species,
                        new_segment_1[i].gene_id,
                        new_segment_2[i].species,
                        new_segment_2[i].gene_id]

                gene.active = False
                gene_family = gene.gene_family

                self.gene_families_second[gene_family].genes.append(new_segment_1[i])
                self.gene_families_second[gene_family].genes.append(new_segment_2[i])
                self.gene_families_second[gene.gene_family].register_event(time, "D", ";".join(map(str, nodes)))

        def make_loss_divisions(time, event):

            lineage = event.lineage
            chromosome = [x for x in self.all_genomes_second[lineage]][0] # This should be corrected, only works if the genomes have a single chromosome
            
            chromosome.obtain_flankings()
            chromosome.obtain_locations()
            
            pseudo = event.pseudo

            c1 = event.sbpL
            c2 = event.sbpR

            r = chromosome.return_affected_region(event.sbpL, event.sbpR, RIGHT)

            if r == None:
                return None

            gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r

            segment = chromosome.obtain_segment(gpositions)
            intergene_segment = chromosome.obtain_intergenic_segment(igpositions[1:])

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

            pseudo_intergenes = []
            if pseudo:
                pseudo_intergenes = intergene_segment

            loss = Loss(int1, int2, c1, c2, specificlen, totallen, lineage, time,
                        pseudo, pseudo_intergenes)
            chromosome.event_history.append(loss)

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
                self.gene_families_second[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))


        def make_inversion_divisions(time, event):

            lineage = event.lineage
            chromosome = [x for x in self.all_genomes_second[lineage]][0] # This should be corrected, only works if the genomes have a single chromosome
            
            chromosome.obtain_flankings()
            chromosome.obtain_locations()

            c1 = event.sbpL
            c2 = event.sbpR

            r = chromosome.return_affected_region(event.sbpL, event.sbpR, RIGHT)

            if r == None:
                return None    

            gpositions, igpositions, leftlengths, rightlengths, int1, int2 = r

                # Get lengths from last intergene before modifying chromosome.
            specificlen = chromosome.intergenes[-1].specific_flanking[1]
            totallen = chromosome.intergenes[-1].total_flanking[1]

            segment = chromosome.obtain_segment(gpositions)
            chromosome.invert_segment(gpositions)

            sleftlen, srightlen, tleftlen, trightlen = 0, 0, 0, 0
            if gpositions[0] > gpositions[-1]:      #The inversion wraps:
                sleftlen, srightlen, tleftlen, trightlen = \
                    chromosome.inversion_wrap_lengths(gpositions)

            inv = Inversion(int1, int2, c1, c2, specificlen, totallen, sleftlen,
                            tleftlen, srightlen, trightlen, lineage, time)
            
            chromosome.invert_divisions(c1,c2) # This function receives the two cuts and inverts all the divisions between them KRISTER

            chromosome.event_history.append(inv)

            scar1 = chromosome.intergenes[igpositions[0]]
            scar2 = chromosome.intergenes[igpositions[-1]]

            scar1.length = leftlengths[0] + rightlengths[0]
            scar2.length = rightlengths[1] + leftlengths[1]

            assert scar1.length == len(inv.afterL)
            assert scar2.length == len(inv.afterR)      

            for gene in segment:
                self.gene_families_second[gene.gene_family].register_event(str(time), "I", ";".join(map(str,[lineage, gene.gene_id])))

        ######

    
        self.all_genomes_second = dict()
        self.gene_families_second = self.initial_gene_families
        self.all_genomes_second["Initial"] = self.initial_genome
        self.all_genomes_second["Root"] = copy.deepcopy(self.initial_genome)

       # We create a list of all the events (T and G) that we will order by time

        all_events = list()
        all_events += [("T",x) for x in self.tree_events]
        for node in self.complete_tree.traverse():
           genome = self.all_genomes[node.name]            
           for chromosome in genome:
               for event in chromosome.event_history: 
                   all_events.append(("G", (event.time, event, chromosome)))
 
       # We start with the initial genome, that we preserved in a copy of the initial genome in the main simulation
       
        genome = self.initial_genome
        
        for items in sorted(all_events, key=lambda x: float(x[1][0])): # This sorts the events by the time

            # We unpack the events

            if items[0] == "T":
                time, etype, lineages = items[1] # In the case that it is a species level event
            else:
                time, event, chromosome = items[1]  # In the case that it is a genome level event
                etype = event.etype
            
            # Species level events

            if etype == "S":
                make_speciation_divisions(time, lineages)
            if etype == "E":
                make_extinction_divisions(time, lineages)
            if etype == "F":
                make_end_divisions(time, lineages)

            # Genome level events

            if etype == "D":
                make_duplication_divisions(time, event)
            if etype == "L":
                make_loss_divisions(time, event)
            if etype == "I":
                print(event)
                c1 = event.sbpL
                c2 = event.sbpR
                print(c1, c2)
                make_inversion_divisions(time, event)


    def write_division_trees(self, division_tree_folder):

        """
        This function writes the division trees and a file documenting the lenghts of every
        division in the initial genome
        
        """

        if not os.path.isdir(division_tree_folder):
            os.mkdir(division_tree_folder)

        for division_family_name, division_family in self.all_division_families.items():
            complete_tree, pruned_tree, _ = division_family.generate_tree()

            with open(os.path.join(division_tree_folder, division_family_name + "_completetree.nwk"), "w") as f:
                f.write(complete_tree)

            if pruned_tree != None:
                with open(os.path.join(division_tree_folder, division_family_name + "_prunedtree.nwk"), "w") as f:
                    f.write(pruned_tree)
            
        with open(os.path.join(division_tree_folder, "Division_lengths.tsv"), "w") as f:
            for division_family_name, division_family in self.all_division_families.items():
                f.write("\t".join(list(map(str,[division_family_name, len(division_family)]))) + "\n")

    





     