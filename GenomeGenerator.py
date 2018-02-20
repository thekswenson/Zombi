import ete3
import os
import AuxiliarFunctions as af
from RatesManager import GeneEvolutionRates
import numpy
import random
from globals import *
from Genome import Genome
import copy
class FamilyOriginator():

    def __init__(self, tree_file, events_file):

        # We need the branch length AND the time of origination of the species

        self.branch_origin = dict()

        with open(events_file) as f:
            f.readline()
            for line in f:
                time, event, node, clades = line.strip().split("\t")
                if event == "SP":
                    c1, c2 = clades.split("+")
                    self.branch_origin[c1] = int(time)
                    self.branch_origin[c2] = int(time)

        with open(tree_file) as f:
            self.whole_tree = ete3.Tree(f.readline().strip(),format=1)

        self.branch_length = dict()

        for node in self.whole_tree.iter_descendants():

            self.branch_length[node.name] = int(node.dist)

        self.vector_names = [x for x in self.branch_length.keys()]
        self.vector_lengths = [self.branch_length[x] for x in self.vector_names]

    def create_families(self):

        # We select first a branch

        branch = numpy.random.choice(self.vector_names, 1, p=af.normalize(self.vector_lengths))[0]

        # We select from an uniform distribution a position for the branch

        time_in_branch = numpy.random.randint(0, self.branch_length[branch])

        # We give the absolute time

        return branch, time_in_branch + self.branch_origin[branch]

SEED = 237
random.seed(SEED)
numpy.random.seed(SEED)

class GeneFamilySimulator():

    def __init__(self, parameters_file, events_file, lineages_in_time_file):

        self.parameters = dict()

        with open(parameters_file) as f:
            for line in f:
                parameter, value = line.strip().split("\t")
                self.parameters[parameter] = value

        self.RM = GeneEvolutionRates(self.parameters)
        self.tree_events = dict()
        self.rates = self.RM.mode_0()

        with open(events_file) as f:
            f.readline()
            for line in f:
                time, dn, ln, cld = line.strip().split("\t")

                if dn == "END":
                    self.total_time = int(time)
                    continue

                elif int(time) not in self.tree_events:
                    self.tree_events[int(time)] = list()
                self.tree_events[int(time)].append((dn,ln,cld))

        self.lineages_in_time = dict()

        with open(lineages_in_time_file) as f:
            f.readline()
            for line in f:
                k,v = line.strip().split("\t")
                self.lineages_in_time[int(k)] = v.split(";")

    def choose_event(self, duplication, transfer, loss):

        draw = numpy.random.choice(["D", "T", "L"], 1, p=af.normalize([duplication, transfer, loss]))
        return draw

    def choose_recipient(self, time_counter, donor, strategy):

        possible_recipients = [x for x in self.lineages_in_time[time_counter] if x != donor]
        if len(possible_recipients) > 1:  # The transfer can go to any other leave, but it cannot be transfer
            recipient = random.choice(possible_recipients)
            return recipient
        else:
            return None

    def evolve_gene_family(self, duplication, transfer, loss, time_counter):

        total_probability_of_event = duplication + transfer + loss

        active_lineages_gt = [x for x in self.gene_family["Gene_tree"].get_leaves() if
                              x.is_alive == True]

        for g_node in active_lineages_gt:

            if g_node.dist == 0:
                continue

            elif numpy.random.uniform(0, 1) <= total_probability_of_event:  # An event takes place

                event = self.choose_event(duplication, transfer, loss)

                if event == "D":
                    self.get_duplicated(g_node, time_counter)

                elif event == "T":

                    recipient = self.choose_recipient(time_counter, g_node.current_branch, 0)
                    if recipient != None:
                        self.get_transferred(g_node, recipient, time_counter)

                elif event == "L":
                    self.get_lost(g_node, time_counter)

    def run_mode_0(self):

        # If rates are global we don't have to compute new each time rates
        # We resort to the global rates already computed

        duplication, transfer, loss = self.rates

        origin_time = self.gene_family["Origin_time"]

        for time_counter in range(origin_time, int(self.total_time)):

            self.increase_distances()

            active_branches = {x.current_branch for x in self.gene_family["Gene_tree"].get_leaves()}

            if time_counter in self.tree_events:

                for event, snode, children in self.tree_events[time_counter]:

                    if snode in active_branches:

                        if event == "EX":
                            self.get_extinct(snode)

                        elif event == "SP":
                            sc1, sc2 = children.split("+")
                            self.gene_tree_speciation(snode, sc1, sc2)

            self.evolve_gene_family(duplication, transfer, loss, time_counter)


    def run_mode_1(self):

        # Rates are family wise. This is control by the rate manager
        # I must keep track of the rates computed with a new type of file

        duplication, transfer, loss = self.RM.mode_0()

        origin_time = self.gene_family["Origin_time"]

        for time_counter in range(origin_time, self.total_time):

            self.increase_distances()

            active_branches = {x.current_branch for x in self.gene_family["Gene_tree"].get_leaves()}

            if time_counter in self.tree_events:

                for event, snode, children in self.tree_events[time_counter]:

                    if snode in active_branches:

                        if event == "EX":
                            self.get_extinct(snode, time_counter)

                        elif event == "SP":
                            sc1, sc2 = children.split("+")
                            self.gene_tree_speciation(snode, sc1, sc2, time_counter)

            self.evolve_gene_family(duplication, transfer, loss, time_counter)

    def run_mode_3(self):

        # Rates are lineage wise. This is control by the rate manager
        pass



    def increase_distances(self):

        active_lineages_gt = [x for x in self.gene_family["Gene_tree"].get_leaves() if
                                  x.is_alive == True]

        for node in active_lineages_gt:
            node.dist += 1

    def origination(self, branch, time_counter, name):

        self.gene_family = dict()
        self.gene_family["Name"] = name
        self.gene_family["Origin_time"] = time_counter
        self.gene_family["Events"] = list()
        self.gene_family["Transfers"] = list()
        self.gene_family["Gene_tree"] = ete3.Tree()
        self.gene_family["Gene_tree"].add_feature("is_alive", True)
        self.gene_family["Gene_tree"].add_feature("current_branch", branch)
        self.gene_family["Gene_tree_extant"] = ete3.Tree()
        self.gene_family["Events"].append(("Origination",time_counter, branch))
        self.gene_family["is_alive"] = True # To avoid too large families
        self.gene_family["Profile"] = dict()

        if time_counter == 0:
            self.gene_family["Profile"]["Root"] = 1

        # To store the events

        self.gene_family["Duplications"] = dict()
        self.gene_family["LeavingTransfers"] = dict()
        self.gene_family["ArrivingTransfers"] = dict()
        self.gene_family["Losses"] = dict()

    def get_duplicated(self, leaf, time_counter):

        self.gene_family["Events"].append(
            ("Duplication", time_counter, leaf.current_branch))

        if leaf.current_branch not in self.gene_family["Duplications"]:
            self.gene_family["Duplications"][leaf.current_branch] = 0

        self.gene_family["Duplications"][leaf.current_branch] +=1

        gc1 = leaf.add_child(dist=0)
        gc1.add_feature("is_alive", True)
        gc1.add_feature("current_branch", leaf.current_branch)

        gc2 = leaf.add_child(dist=0)
        gc2.add_feature("is_alive", True)
        gc2.add_feature("current_branch", leaf.current_branch)

        leaf.name = "D"
        leaf.is_alive = False

    def get_lost(self, leaf, time_counter):

        self.gene_family["Events"].append(
            ("Loss", time_counter, leaf.current_branch))

        if leaf.current_branch not in self.gene_family["Losses"]:
            self.gene_family["Losses"][leaf.current_branch] = 0

        self.gene_family["Losses"][leaf.current_branch] +=1

        leaf.is_alive = False
        leaf.name = "L"

    def get_transferred(self, leaf, recipient, time_counter):

        rt = float(self.parameters["REPLACEMENT_T"])

        self.gene_family["Events"].append(
            ("LeavingTransfer", time_counter, leaf.current_branch))
        self.gene_family["Events"].append(
            ("ArrivingTransfer", time_counter, recipient))

        self.gene_family["Transfers"].append((leaf.current_branch, recipient))

        if leaf.current_branch not in self.gene_family["LeavingTransfers"]:
            self.gene_family["LeavingTransfers"][leaf.current_branch] = 0

        if leaf.current_branch not in self.gene_family["ArrivingTransfers"]:
            self.gene_family["ArrivingTransfers"][leaf.current_branch] = 0

        self.gene_family["LeavingTransfers"][leaf.current_branch] += 1
        self.gene_family["ArrivingTransfers"][leaf.current_branch] += 1

        leaf.name = "TRANSFER"

        gc1 = leaf.add_child(dist=0)
        gc1.add_feature("is_alive", True)
        gc1.add_feature("current_branch", recipient)
        gc1.name = recipient

        gc2 = leaf.add_child(dist=0)
        gc2.add_feature("is_alive", True)
        gc2.add_feature("current_branch", leaf.current_branch)
        gc2.name = leaf.current_branch

        leaf.is_alive = False

        other_copies = [x for x in self.gene_family["Gene_tree"].get_leaves() if x.current_branch == recipient and x != 0 ]

        if len(other_copies) != 0 and numpy.random.uniform(0,1) < rt:

            # So far, this is like a normal transfer, but now I have to kill one of the extant genes in that branch
            replaced_lineage = random.choice(other_copies)
            replaced_lineage.is_alive = False

            # I need to add the gene lost by replacement

            self.gene_family["Events"].append(
                ("LossByReplacement", time_counter, recipient))

    def get_extinct(self, sp_branch):

        g_leaves = self.gene_family["Gene_tree"].get_leaves()

        #self.gene_family["Events"].append(
        #    ("Extinction", time_counter * TIME_INCREASE, sp_branch))

        self.gene_family["Profile"][sp_branch] = 0

        for g_leaf in g_leaves:

            if g_leaf.current_branch == sp_branch and g_leaf.is_alive == True:

                self.gene_family["Profile"][sp_branch] += 1

                g_leaf.is_alive = False

    def gene_tree_speciation(self, sp, c1, c2):

        #self.gene_family["Events"].append(
        #    ("Speciation", time_counter * TIME_INCREASE, sp + "->" + c1 + ";" + c2))

        g_leaves = self.gene_family["Gene_tree"].get_leaves()

        self.gene_family["Profile"][sp] = 0

        for g_leaf in g_leaves:

            if g_leaf.current_branch == sp and g_leaf.is_alive == True:

                self.gene_family["Profile"][sp] += 1

                gc1 = g_leaf.add_child(dist=0)
                gc1.add_feature("is_alive", True)
                gc1.add_feature("current_branch", c1)

                gc2 = g_leaf.add_child(dist=0)
                gc2.add_feature("is_alive", True)
                gc2.add_feature("current_branch", c2)

                g_leaf.name = "SP" + sp
                g_leaf.is_alive = False

    def complete_gene_family_information(self):

        # Add root length

        myroot = self.gene_family["Gene_tree"].get_tree_root()
        one_leaf, distance_to_present = myroot.get_farthest_leaf()
        myroot.dist = (self.total_time - distance_to_present) - self.gene_family["Origin_time"]

        # Add names

        genetree = self.gene_family["Gene_tree"]
        names = dict()

        for leaf in genetree.get_leaves():

            leaf.name = leaf.current_branch
            myname = leaf.name.split("_")[0]

            if myname in names:
                names[myname] += 1
            else:
                names[myname] = 1

            if leaf.is_alive:
                leaf.name = myname + "_" + "A" + "_" + str(names[leaf.name])
            else:
                leaf.name = myname + "_" + "E" + "_" + str(names[leaf.name])

        for leaf in genetree.get_leaves():
            if not leaf.is_alive:
                continue
            if leaf.current_branch not in self.gene_family["Profile"]:
                self.gene_family["Profile"][leaf.current_branch] = 0
            self.gene_family["Profile"][leaf.current_branch] += 1

    def get_gene_family_tree(self):

        return self.gene_family["Gene_tree"].write(format=1)

    def output_profile(self, profile, which):

        with open(profile, "r+") as f:
            header = f.readline().strip().split("\t")[1:]
            myline = list()
            myline.append(str(self.gene_family["Name"]))

            for name in header:

                try: myline.append(str(self.gene_family[which][name]))
                except: myline.append(str(0))

            myline = "\t".join(myline) +"\n"
            f.write(myline)

    def write_transfers(self, transfers_file):

        with open(transfers_file, "a") as f:
            for dn, rc in self.gene_family["Transfers"]:
                line = self.gene_family["Name"] + "\t" + dn + "\t" + rc + "\n"
                f.write(line)

    def write_log(self, events_logfile):

        with open(events_logfile, "w") as f:

            f.write("Time\tEvent\tNode\n")

            for event, time, node in self.gene_family["Events"]:

                myline = "\t".join(map(str,[event,time,node])) + "\n"
                f.write(myline)



class GenomeSimulator():

    def __init__(self, parameters_file, events_file, lineages_in_time_file):

        self.homologous = dict()

        self.parameters = dict()
        self.species_tree = ete3.Tree()

        with open(parameters_file) as f:
            for line in f:
                parameter, value = line.strip().split("\t")
                self.parameters[parameter] = value

        self.tree_events = dict()
        self._start_tree()

        with open(events_file) as f:
            f.readline()
            for line in f:
                time, dn, ln, cld = line.strip().split("\t")
                if dn == "END":
                    self.total_time = int(time)
                    continue
                elif int(time) not in self.tree_events:
                    self.tree_events[int(time)] = list()
                self.tree_events[int(time)].append((dn,ln,cld))

        self.lineages_in_time = dict()

        with open(lineages_in_time_file) as f:
            f.readline()
            for line in f:
                k,v = line.strip().split("\t")
                self.lineages_in_time[int(k)] = v.split(";")

    def _start_tree(self):

        gnm = Genome()
        gnm.start_genome(int(self.parameters["STEM_FAMILIES"]))
        root = self.species_tree.get_tree_root()
        root.name = "Root"
        root.add_feature("Genome",gnm)
        root.add_feature("Is_alive", True)

        for i in range(1, int(self.parameters["STEM_FAMILIES"])+1):
            self.homologous[str(i)] = dict()
            self.homologous[str(i)]["Copies"] = 1
            self.homologous[str(i)]["Events"] = list()

    def choose_event(self, duplication, transfer, loss, inversion, translocation, origination):

        draw = numpy.random.choice(["D", "T", "L", "I", "C", "O"], 1, p=af.normalize([duplication, transfer, loss, inversion, translocation, origination]))
        return draw

    def choose_recipient(self, time_counter, donor, strategy):

        possible_recipients = [x for x in self.lineages_in_time[time_counter] if x != donor]
        if len(possible_recipients) > 1:
            recipient = random.choice(possible_recipients)
            return recipient
        else:
            return None

    def evolve_genomes(self, duplication, transfer, loss, inversion, translocation, origination, time_counter):

        total_probability_of_event = duplication + transfer + loss + inversion + translocation

        active_genomes = [x for x in self.species_tree.get_leaves() if x.is_alive == True]
        random.shuffle(active_genomes)

        for node in active_genomes:

            genome = node.Genome

            if numpy.random.uniform(0, 1) <= total_probability_of_event:  # An event takes place

                event = self.choose_event(duplication, transfer, loss, inversion, translocation, origination)

                if event == "D":
                    a = genome.obtain_affected_genes()
                    genome.duplicate_segment(self.homologous, time_counter, a)

                elif event == "T":
                    recipient = self.choose_recipient(time_counter, node.name, 0)

                    if recipient == None:
                        continue

                    a = genome.obtain_affected_genes()
                    segment = genome.obtain_segment(a)

                    old_segment = list()
                    new_segment = list()

                    for i in a:
                        sense, gf, cp = genome.genes[i].split("_")
                        self.homologous[gf]["Copies"] += 1
                        name1  = sense + "_" + gf + "_" + str(self.homologous[gf]["Copies"])
                        old_segment.append(name1)
                        self.homologous[gf]["Copies"] += 1
                        name2 = sense + "_" + gf + "_" + str(self.homologous[gf]["Copies"])
                        new_segment.append(name2)
                        self.homologous[gf]["Events"].append(
                            ("T", str(time_counter), cp + "_" + name1.split("_")[2] + "_" + name2.split("_")[2]))

                    # Now I have prepared the two segments. First I am going to update de donor segment

                    elements_to_remove = [genome.genes[x] for x in a]

                    for element in elements_to_remove:
                        genome.genes.remove(element)

                    position = a[0]

                    for i, x in enumerate(old_segment):
                        genome.genes.insert(position + i + 1, x)

                    # Then I update the receptor segment

                        my_recipient = self.species_tree&recipient
                        recipient_genome = my_recipient.Genome
                        p = recipient_genome.select_random_position()
                        recipient_genome.insert_segment(p, new_segment)

                elif event == "L":

                    # We have to check that the minimal size has not been attained

                    if len(genome.genes) <= int(self.parameters["MIN_GENOME_SIZE"]):
                        continue

                    a = genome.obtain_affected_genes()
                    genome.loss_segment(self.homologous, time_counter, a)

                elif event == "I":
                    a = genome.obtain_affected_genes()
                    genome.invert_segment(self.homologous,time_counter, a)

                elif event == "C":
                    a = genome.obtain_affected_genes()
                    genome.translocate_segment(self.homologous, time_counter, a)

                elif event == "O":
                    a = genome.obtain_affected_genes()
                    genome.invert_segment(self.homologous, time_counter, a)

    def run(self):

        duplication, transfer, loss, inversion, translocation, origination = (0.001, 0.000, 0.0001, 0.00, 0.000, 0.00)

        for time_counter in range(int(self.total_time + 1)):

            if time_counter in self.tree_events:

                for event, snode, children in self.tree_events[time_counter]:

                    if event == "EX":
                        self.get_extinct(time_counter, snode)

                    elif event == "SP":

                        # First we write the ancestral genome

                        snode.Genome.write_genome()

                        sc1, sc2 = children.split("+")
                        self.get_speciated(time_counter, snode, sc1, sc2)

            self.evolve_genomes(duplication, transfer, loss, inversion, translocation, origination, time_counter)
            self.increase_distances()

        print(time_counter)
        print(self.species_tree.write(format=1))

        survivors = [x for x in self.species_tree.get_leaves() if x.is_alive == True]
        for n in survivors:
            print(n.Genome.genes)

    def increase_distances(self):

        active_lineages = [x for x in self.species_tree.get_leaves() if x.is_alive == True]

        for node in active_lineages:
            node.dist += 1


    def get_extinct(self, time, sp):

        sp = self.species_tree&sp
        sp.is_alive = False
        parent_genome = sp.Genome

        for i, gene in enumerate(parent_genome.genes):
            sense, gf, hml = gene.split("_")
            self.homologous[gf]["Events"].append(("E", time, hml))

    def get_speciated(self, time, sp, c1, c2):

        sp = self.species_tree&sp
        parent_genome = sp.Genome

        sc1 = sp.add_child(dist=0)
        sc1.name = c1
        sc1.add_feature("is_alive", True)
        sc1.add_feature("Genome", copy.deepcopy(parent_genome))
        genes_affected_1 = sc1.Genome.update_homologous(self.homologous)

        sc2 = sp.add_child(dist=0)
        sc2.name = c2
        sc2.add_feature("is_alive", True)
        sc2.add_feature("Genome", copy.deepcopy(parent_genome))
        genes_affected_2 = sc2.Genome.update_homologous(self.homologous)

        for i, gene in enumerate(parent_genome.genes):

            sense, gf, hml = gene.split("_")
            speciation_event = "_".join((sp.name, hml, genes_affected_1[i].split("_")[2], genes_affected_2[i].split("_")[2]))
            self.homologous[gf]["Events"].append(("S", time, speciation_event))

        sp.is_alive = False

    def complete_gene_family_information(self):

        # Add root length

        myroot = self.gene_family["Gene_tree"].get_tree_root()
        one_leaf, distance_to_present = myroot.get_farthest_leaf()
        myroot.dist = (self.total_time - distance_to_present) - self.gene_family["Origin_time"]

        # Add names

        genetree = self.gene_family["Gene_tree"]
        names = dict()

        for leaf in genetree.get_leaves():

            leaf.name = leaf.current_branch
            myname = leaf.name.split("_")[0]

            if myname in names:
                names[myname] += 1
            else:
                names[myname] = 1

            if leaf.is_alive:
                leaf.name = myname + "_" + "A" + "_" + str(names[leaf.name])
            else:
                leaf.name = myname + "_" + "E" + "_" + str(names[leaf.name])

        for leaf in genetree.get_leaves():
            if not leaf.is_alive:
                continue
            if leaf.current_branch not in self.gene_family["Profile"]:
                self.gene_family["Profile"][leaf.current_branch] = 0
            self.gene_family["Profile"][leaf.current_branch] += 1

    def get_gene_family_tree(self):

        return self.gene_family["Gene_tree"].write(format=1)

    def output_profile(self, profile, which):

        with open(profile, "r+") as f:
            header = f.readline().strip().split("\t")[1:]
            myline = list()
            myline.append(str(self.gene_family["Name"]))

            for name in header:

                try: myline.append(str(self.gene_family[which][name]))
                except: myline.append(str(0))

            myline = "\t".join(myline) +"\n"
            f.write(myline)

    def write_transfers(self, transfers_file):

        with open(transfers_file, "a") as f:
            for dn, rc in self.gene_family["Transfers"]:
                line = self.gene_family["Name"] + "\t" + dn + "\t" + rc + "\n"
                f.write(line)

    def write_log(self, events_logfolder):

        for gf in self.homologous:

            with open(os.path.join(events_logfolder, self.parameters["PREFIX"] + str(gf)) + "_events.tsv", "w") as f:

                f.write("Time\tEvent\tNode\n")

                for time,event,node in self.homologous[gf]["Events"]:

                    line = "\t".join(map(str,[time, event, node])) + "\n"
                    f.write(line)

    def generate_gene_tree(self):

        path = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/RESTART/TEST2/ParallelGeneFamilies/FAM2_events.tsv"

        events = dict()

        with open(path) as f:
            f.readline()

            for line in f:
                event, time, nodes = line.strip().split("\t")
                time = int(time)
                if time not in events:
                    events[time] = list()
                events[time].append((event, nodes))

        gene_tree = ete3.Tree()
        root = gene_tree.get_tree_root()
        root.name = "1"
        root.add_feature("is_alive", True)
        root.add_feature("current_branch", "")

        for time in range(self.total_time):

            if time in events:

                for event, nodes in events[time]:

                    if event == "S":

                        sp,parent,c1,c2 = nodes.split("_")

                        n = gene_tree&parent
                        n.is_alive = False
                        n.current_branch = sp

                        gc1 = n.add_child()
                        gc1.name = c1
                        gc1.add_feature("is_alive", True)

                        gc2 = n.add_child()
                        gc2.name = c2
                        gc2.add_feature("is_alive", True)

                    elif event == "E":

                        n = gene_tree&nodes
                        n.is_alive = False

                    elif event == "D":

                        parent, c1, c2 = nodes.split("_")

                        n = gene_tree & parent
                        n.is_alive = False

                        gc1 = n.add_child()
                        gc1.name = c1
                        gc1.add_feature("is_alive", True)

                        gc2 = n.add_child()
                        gc2.name = c2
                        gc2.add_feature("is_alive", True)

                    elif event == "T":

                        parent, c1, c2 = nodes.split("_")

                        n = gene_tree & parent
                        n.is_alive = False

                        gc1 = n.add_child()
                        gc1.name = c1
                        gc1.add_feature("is_alive", True)

                        gc2 = n.add_child()
                        gc2.name = c2
                        gc2.add_feature("is_alive", True)

                    elif event == "L":

                        n = gene_tree&nodes
                        n.is_alive = False



            active_lineages = [x for x in gene_tree.get_leaves() if x.is_alive == True]

            for node in active_lineages:
                node.dist+=1

        print(gene_tree.write(format=1))










































