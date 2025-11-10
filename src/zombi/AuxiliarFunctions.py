import glob
import re
import ete3
import numpy
import sys
import scipy
import scipy.stats as ss
import pandas as pd

from collections import defaultdict
from itertools import pairwise
from typing import Any
from pathlib import Path
from numpy.random import Generator as npGenerator
from BCBio import GFF
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from .Events import FER, LOSS, ORIG, TDUP
from .Filenames import COMPLETEsuffix, GENOMEsuffix, SUBSTITUTIONSCALEDsuffix
from .Filenames import PRUNEDsuffix, SAMPLEDsuffix, INITIALGENOMEINFO
from .Filenames import GENEFAMILYINFO


def normalize(array):
    total = numpy.sum(array)
    if total != 0:
        return (array/total)
    else:
        return array
    
    
def normalize_middle(array):
    middle = array[int(len(array)/2)]
    
    if middle != 0:
        return array/middle
    else:
        return array
    

def inverse(array):
    transformed_array = 1./numpy.array(array)
    return (transformed_array)


def logtransform(array):
    transformed_array = numpy.log(array)
    return (transformed_array)

def calculate_mean_genome_size(myfile):
    genomes = list()
    with open(myfile) as f:
        header = f.readline().strip().split("\t")[1:]
        for node in header:
            genomes.append(0)
        for line in f:
            handle = line.strip().split("\t")[1:]
            for i,gf in enumerate(handle):
                genomes[i] += int(gf)
    return(numpy.array(genomes).mean())

def divide_by_time_increase(array):
    transformed_array = numpy.log(array)
    return (transformed_array)

def read_parameters(parameters_file) -> dict[str, str]:

    parameters = dict()

    with open(parameters_file) as f:

        try:
            for i, line in enumerate(f):

                if line[0] == "#" or line == "\n":
                    continue

                if "\t" in line and " " in line:
                    sys.exit(f'ERROR: parameter file "{parameters_file}" has a '
                             f'line\n"{line.strip()}" containing a mix of tabs '
                             f'and spaces!')
                elif "\t" in line:
                    parameter, value = line.strip().split("\t")
                    parameters[parameter] = value
                elif " " in line:
                    parameter, value = line.strip().split(" ")
                    parameters[parameter] = value

        except Exception:
            sys.exit(f'ERROR: problem reading line {i+1} of '
                     f'"{parameters_file}":\n"{line.rstrip()}"')

    return parameters


def read_seed(parameters_file: Path):

    myseed = 0

    with open(parameters_file) as f:

        for line in f:

            if line[0] == "#" or line == "\n":
                continue
            if "SEED" not in line:
                continue

            if "\t" in line:
                _, value = line.strip().split("\t")
                myseed = int(value)
            elif " " in line:
                _, value = line.strip().split(" ")
                myseed = int(value)

    return myseed

def read_empirical_rates(rates_file: str, scale_rates = 1.0):

    empirical_rates = list()

    with open(rates_file) as f:
        f.readline()
        for line in f:
            _, d, t, l = line.strip().split("\t")
            d, t, l = [x/scale_rates for x in map(float,[d, t, l])]
            empirical_rates.append((float(d), float(t), float(l)))

    return empirical_rates


def obtain_value(value, nprngen: npGenerator) -> float:

    handle = value.split(":")

    if handle[0] == "f":
        # Fixed value
        value =  float(handle[1])

    elif handle[0] == "n":
        # normal distribution
        params = handle[1].split(";")
        value = abs(nprngen.normal(float(params[0]), float(params[1])))

    elif handle[0] == "l":
        # lognormal distribution
        params = handle[1].split(";")
        value = abs(nprngen.lognormal(float(params[0]), float(params[1])))

    elif handle[0] == "u":
        # uniform distribution
        params = handle[1].split(";")
        value = abs(nprngen.uniform(float(params[0]), float(params[1])))

    elif handle[0] == "g":
        # geometric distribution
        value = float(nprngen.geometric(float(handle[1])))

    elif handle[0] == "e":
        # exponential distribution
        value = nprngen.exponential(float(handle[1]))

    return value

def discretize(alpha, ncat, disttype="lognorm"):

    # adapted from Kevin Gori
    # Taken from https://gist.github.com/kgori/95f604131ce92ec15f4338635a86dfb9

    if disttype == "gamma":
        dist = ss.gamma(alpha, scale=1 / alpha)
    elif disttype == "lognorm":
        dist = ss.lognorm(s=alpha, scale=numpy.exp(0.5 * alpha**2))
   
    quantiles = dist.ppf(numpy.arange(0, ncat) / ncat)    
    rates = numpy.zeros(ncat, dtype=numpy.double)

    assert isinstance(dist, ss.rv_continuous)
    for i in range(ncat-1):
        rates[i] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), 
                                               quantiles[i], quantiles[i+1])[0]
    rates[ncat-1] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), 
                                                quantiles[ncat-1], numpy.inf)[0]
    return rates

def sample_from_dirichlet(n, nprngen: npGenerator):
    return nprngen.dirichlet([1] * n)


def prepare_sequence_parameters(parameters_file: Path) -> dict[str, Any]:
    """
    Convert some of the S parameters to the correct type.
    """
    parameters: dict[str, Any] = read_parameters(parameters_file)
    for parameter, value in parameters.items():

        if parameter == "SEED":
            parameters[parameter] = int(value)
            if parameters[parameter] == 0:
                parameters[parameter] = None

        if parameter == "SEQUENCE_SIZE" or parameter == "VERBOSE" or parameter == "SEED":
            parameters[parameter] = int(value)

        if parameter == "SCALING":
            parameters[parameter] = float(value)

    return parameters


def prepare_species_tree_parameters(parameters: dict[str, Any],
                                    nprngen: npGenerator) -> dict[str, Any]:
    """
    Convert some of the T parameters to the correct type.
    """
    for parameter, value in parameters.items():

        if parameter == "TURNOVER":
            parameters[parameter] = obtain_value(value, nprngen)

        elif parameter == "TOTAL_TIME" or parameter == "SCALE_TREE":
            parameters[parameter] = float(value)

        elif parameter == "LINEAGE_PROFILE":
            parameters[parameter] = [tuple([int(j) for j in x.split("-")]) for x in value.split(";")]
            
        elif parameter == "MASSIVE_EXTINCTION":
            parameters[parameter] = [tuple([float(j) for j in x.split("-")]) for x in value.split(";")]
            
        elif(parameter == "SPECIES_EVOLUTION_MODE" or parameter == "N_LINEAGES" or parameter == "MIN_LINEAGES" or
             parameter == "TOTAL_LINEAGES" or parameter == "STOPPING_RULE" or parameter == "MAX_LINEAGES" or
             parameter == "VERBOSE" or parameter == "SEED" or
             parameter == "NUM_SPECIATION_RATE_CATEGORIES" or parameter == "NUM_EXTINCTION_RATE_CATEGORIES" or
             parameter == "SIMULATE_SEQUENCES" or parameter == "SCALE_GENE_TREES"):
            parameters[parameter] = int(value)

        elif(parameter == "SPECIATION" or
             parameter == "EXTINCTION" or
             parameter == "MASS_EXTINCTION" or
             parameter == "SHIFT_SPECIATION_RATE_FREQUENCY" or
             parameter == "SHIFT_EXTINCTION_RATE_FREQUENCY" or
             parameter == "BASE_EXTINCTION" or
             parameter == "BASE_SPECIATION"):
            pass

        else:
            sys.exit(f'ERROR: unknown parameter "{parameter}" given to tree '
                     f'simulator.\n       Did you give the correct file?')

    return parameters


def get_complementary_sequence(sequence: str) -> str:

    new_sequence = sequence.replace("A","x").replace("T","y").replace("C","v").replace("G","w")
    new_sequence = new_sequence.replace("x","T").replace("y","A").replace("v","G").replace("w","C")
    new_sequence = new_sequence[::-1]
    return new_sequence


def fasta_reader(fasta_file: Path):

    with open(fasta_file) as f:

        seq = ""
        for line in f:
            if ">" == line[0]:
                if seq != "":
                    yield header, seq
                    header = line.strip()
                    seq = ""
                else:
                    header = line.strip()
                    seq = ""
            else:
                seq += line.strip()

        yield header, seq

def fasta_writer(outfile: Path, entries):

    x = 80
    with open(outfile, "w") as f:
        for h, seq in entries:
            f.write(h + "\n")
            lines = [seq[i: i + x] for i in range(0, len(seq), x)]
            for line in lines:
                f.write(line +"\n")


def prepare_genome_parameters(parameters_file: Path):

    parameters: dict[str, Any] = read_parameters(parameters_file)
    for parameter, value in parameters.items():

        #if parameter == "DUPLICATION_EXTENSION" or parameter == "TRANSFER_EXTENSION" \
        #        or parameter == "LOSS_EXTENSION" or parameter == "INVERSION_EXTENSION" or \
        #        parameter == "TRANSPOSITION_EXTENSION" or parameter == "ORIGINATION_EXTENSION":

        #    parameters[parameter] = obtain_value(value)


        if parameter == "ROOT_GENOME":
            parameters[parameter] = value.split(";")

        if(parameter == "REPLACEMENT_TRANSFER" or parameter == "ALPHA" or
           parameter == "SCALE_TREE"):
            parameters[parameter] = float(value)

        if(parameter == "PROFILES" or parameter == "EVENTS_PER_BRANCH" or
           parameter == "GENE_TREES" or parameter == "PRUNE_TREES" or
           parameter == "TRANSFER_PREFERENCE" or
           parameter == "RECONCILED_TREES" or parameter == "VERBOSE" or
           parameter == "MIN_GENOME_SIZE" or parameter == "ALL_GENOMES" or
           parameter == "EXTENSION_MULTIPLIER" or
           parameter == "SEED" or parameter == "GENEORDER_EVENTS_PER_BRANCH"):

            parameters[parameter] = int(value)

        if parameter == "ASSORTATIVE_TRANSFER" or parameter == "SCALE_RATES":
            if not (parameters[parameter] == "True" or parameters[parameter] == "False"):
                sys.exit(f'ERROR: parameter "{parameter}" must be either "True" or "False".')
            parameters[parameter] = parameters[parameter] == "True"

    return parameters


def generate_events(tree_file: Path):

    events = []

    with open(tree_file) as f:
        treeline = f.readline().strip()
        tree = ete3.Tree(treeline, format=1)

    root: ete3.TreeNode = tree.get_tree_root()
    root.name = "Root"

    ## There is probably a better way to write this


    try:
        root.dist = float(treeline.split(")")[-1].split(":")[1].replace(";",""))

    except:
        pass

    total_time = root.get_farthest_leaf()[1]
    assert isinstance(total_time, float)

    # There might be slighlty variations in the branch length that we have to account for. So all nodes
    # that are at 0.1% distance of the fathes leaf will be considered to be alive

    error_margin = total_time * 0.001

    nodes = list()

    for node in tree.traverse():   #type: ignore    

        node_dist = node.get_distance(root)

        if node.is_leaf():
            if  total_time <= node_dist + error_margin  and total_time >= node_dist - error_margin:
                nodes.append((node, "A", node_dist))
            else:
                nodes.append((node, "E", node_dist))
        else:
            nodes.append((node, "S", node_dist))

    # We order the events # This can be more efficient

    nodes = sorted(nodes, key = lambda x: x[2])

    for node, estate, time in nodes:
        if estate == "A":
            events.append((str(time + root.dist), "F", node.name))
        elif estate == "E":
            events.append((str(time + root.dist), "E", node.name))
        elif estate == "S":
            c1, c2 = node.get_children()
            events.append((str(time + root.dist), "S", ";".join((node.name, c1.name, c2.name))))

    return events


def return_vector_of_distances(self, tree_file):

    self.distances_to_root = dict()

    with open(tree_file) as f:

        self.mytree = ete3.Tree(f.readline().strip(), format=1)
        root = self.mytree.get_tree_root()
        root.name = "Root"
        for node in self.mytree.traverse():     #type: ignore
            if node.is_root():
                continue

            self.distances_to_root[node.name] = (node, node.get_distance(root))

def choose_advanced_recipient(self, time, alive_lineages, donor,
                              nprngen: npGenerator):

    # Chooses and advanced recipient according to the logarithm of the phylogenetic distance

    possible_recipients = list()
    weights = list()

    mydonor = self.distances_to_root[donor][0]

    for recipient in alive_lineages:

        if donor == recipient:
            continue

        myrecipient = self.distances_to_root[recipient][0]
        phylo_d = mydonor.get_distance(myrecipient)

        td = phylo_d + (2 * time) - self.distances_to_root[donor][1] - self.distances_to_root[recipient][1]

        possible_recipients.append(recipient)
        weights.append(td)

    draw = nprngen.choice(possible_recipients, 1, p= normalize(weights))
    raise NotImplementedError("This appears to be unfinished")


def generate_newick_trees(events):

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

    surviving_nodes = dict()
    times = dict()

    for current_time, event, nodes in reversed(events):

        if event == "F":

            times[nodes] = float(current_time)
            surviving_nodes[nodes] = {"state": 1, "descendant": "None"}

        elif event == "E":

            times[nodes] = float(current_time)
            surviving_nodes[nodes] = {"state": 0, "descendant": "None"}

        elif event == "S":

            p, c1, c2 = nodes.split(";")

            times[p] = float(current_time)

            if surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 1:

                surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + c2}

            elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 0:

                surviving_nodes[p] = {"state": 0, "descendant": "None"}

            elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == -1:

                mynode1 = find_descendant(surviving_nodes, c1)
                mynode2 = find_descendant(surviving_nodes, c2)

                surviving_nodes[p] = {"state": 1, "descendant": mynode1 + ";" + mynode2}


            elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 0:

                surviving_nodes[p] = {"state": -1, "descendant": c1}

            elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 1:

                surviving_nodes[p] = {"state": -1, "descendant": c2}


            elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == -1:

                mynode = find_descendant(surviving_nodes, c2)
                surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + mynode}

            elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 1:

                mynode = find_descendant(surviving_nodes, c1)
                surviving_nodes[p] = {"state": 1, "descendant": mynode + ";" + c2}


            elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 0:

                mynode = find_descendant(surviving_nodes, c1)
                surviving_nodes[p] = {"state": -1, "descendant": mynode}

            elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == -1:

                mynode = find_descendant(surviving_nodes, c2)
                surviving_nodes[p] = {"state": -1, "descendant": mynode}

    extanttree = ete3.Tree()
    completetree = ete3.Tree()
    eroot = extanttree.get_tree_root()
    eroot.name = ""
    wroot = completetree.get_tree_root()
    wroot.name = "Root"

    wquick_nodes = dict()
    equick_nodes = dict()

    wquick_nodes["Root"] = wroot

    for i, values in enumerate(events):

        current_time, event, nodes = values

        if event == "S":

            p, c1, c2 = nodes.split(";")

            mynode = wquick_nodes[p]
            myc1 = mynode.add_child()
            myc2 = mynode.add_child()
            myc1.name = c1
            myc2.name = c2
            myc1.dist = times[c1] - times[p]
            myc2.dist = times[c2] - times[p]

            wquick_nodes[c1] = myc1
            wquick_nodes[c2] = myc2

            state = surviving_nodes[p]["state"]

            if state == 1:  # Now the extant tree

                c1name, c2name = surviving_nodes[p]["descendant"].split(";")

                if eroot.name == "":
                    eroot.name = p
                    equick_nodes[p] = eroot

                mynode = equick_nodes[p]

                myc1 = mynode.add_child()
                myc2 = mynode.add_child()

                myc1.name = c1name
                myc2.name = c2name

                myc1.dist = times[c1name] - times[p]
                myc2.dist = times[c2name] - times[p]

                equick_nodes[c1name] = myc1
                equick_nodes[c2name] = myc2

    return completetree.write(format=1, format_root_node=True), extanttree.write(format=1, format_root_node=True)



def generate_gene_tree(events):

 # THIS IS FOR GENERATING GENE TREES


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

    surviving_nodes = dict()
    times = dict()

    for current_time, event, nodes in reversed(events):

        if event == "F":

            nodename = nodes.replace(";","_")

            times[nodename] = float(current_time)
            surviving_nodes[nodename] = {"state": 1, "descendant": "None"}

        elif event == "E" or event == LOSS:

            nodename = nodes.replace(";", "_")

            times[nodename] = float(current_time)
            surviving_nodes[nodename] = {"state": 0, "descendant": "None"}

        elif event == "S" or event == TDUP or event == FER:

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



    extanttree = ete3.Tree()
    completetree = ete3.Tree()
    eroot = extanttree.get_tree_root()
    eroot.name = ""

    wquick_nodes = dict()
    equick_nodes = dict()

    for i, values in enumerate(events):

        current_time, event, nodes = values

        if event == ORIG:

            wroot = completetree.get_tree_root()
            wroot.name = nodes + "_1"
            wquick_nodes[wroot.name] = wroot

        if event == "S" or event == TDUP or event == FER:

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


    if len(completetree) == 0:
        completetree = ";"

    elif len(completetree) == 1:
        completetree =  completetree.get_leaves()[0].name + ";"

    else:
        completetree = completetree.write(format=1, format_root_node=True)

    if len(extanttree) == 0:
        extanttree = ";"

    elif len(extanttree) == 1:
        extanttree = extanttree.get_leaves()[0].name + ";"

    else:
        extanttree = extanttree.write(format=1, format_root_node=True)

    return completetree, extanttree


def write_pruned_sequences(tree_file: str, fasta_folder: Path, scaled=False):
    with open(tree_file) as f:
        line = f.readline().strip()
        if "(" not in line or line == ";":
            return None
        else:
            my_tree = ete3.Tree(line, format=1)

    surviving_nodes = {x.name for x in my_tree.get_leaves()}
    tree_num = tree_file.split("/")[-1].split("_")[0]

    if not scaled:
        entries = fasta_reader(fasta_folder / f"{tree_num}{COMPLETEsuffix}")
    else:
        entries = fasta_reader(fasta_folder / f"{tree_num}{SUBSTITUTIONSCALEDsuffix}")

    clean_entries = list()
    for h, seq in entries:
        if h[1:] in surviving_nodes:
            clean_entries.append((h, seq))

    fasta_writer(fasta_folder / f"{tree_num}{PRUNEDsuffix}", clean_entries)


def write_sampled_sequences(tree_file: str, infasta_folder: Path, outfasta_folder: Path):

    with open(tree_file) as f:
        line = f.readline().strip()
        if "(" not in line or line == ";":
            return None
        else:
            my_tree = ete3.Tree(line, format=1)
    surviving_nodes = {x.name for x in my_tree.get_leaves()}
    file_name = tree_file.split("/")[-1].split("_")[0]
    entries = fasta_reader(infasta_folder / f"{file_name}{COMPLETEsuffix}")

    clean_entries = list()
    for h, seq in entries:
        if h[1:] in surviving_nodes:
            clean_entries.append((h, seq))

    fasta_writer(outfasta_folder / f"{file_name}{SAMPLEDsuffix}", clean_entries)


def parse_GFF(gff_file: Path, sort=True) -> tuple[int, list[SeqFeature]]:
    """
    Extract the chromosome length and the genes from the given GFF file.
    Genes are in Biopython SeqFeature format:

        https://biopython.org/docs/latest/api/Bio.SeqFeature.html

    Parameters
    ----------
    gff_file : Path
        the GFF file to parse
    sort : bool, default true
        return the genes sorted by start index

    Returns
    -------
    tuple[int, list[SeqFeature]]
        the chromosome length along with Biopython SeqFeatures for all of the
        genes. The indices for seqfeatures are pythonic (zero indexed, end not
        inclusive).
    """
            #Get the CDSs and genome length:
    #import pprint; pprint.pprint(examiner.available_limits(gff_file))
    genome_len = 0
    genes: list[SeqFeature] = []
    try:
        for rec in GFF.parse(gff_file, target_lines=1000,
                             limit_info={'gff_type': ['CDS']}):
            if 'sequence-region' in rec.annotations:
                genome_len = rec.annotations['sequence-region'][0][2]

            for feature in rec.features:
                if feature.type == 'CDS':
                    genes.append(feature)
    except TypeError as e:
        sys.exit(f'Problem opening GFF file (remove fasta lines?):\n{e}')
    except AttributeError as e:
        sys.exit(f'Problem opening GFF file (remove fasta lines?):\n{e}')

        #Assure that the genes are within the sequence length:
    for gene in genes:
        if gene.location.end > genome_len:                  #type: ignore
            sys.exit(f'There is a problem in the GFF file!  Either\n'
                     f'  1. the sequence length {genome_len} is incorrect, or\n'
                     f'  2. the feature {gene.id} is larger than the length.\n'
                     f'\nWe expect the first line starting with '
                     f'"##sequence-region" to contain the\n'
                     f'genome length (e.g. "##sequence-region L_1 1 4412837").')

    if sort:                    #sort by start index
        genes.sort(key=lambda f: f.location.start)          #type: ignore   

    if not genome_len:
        sys.exit(f'No sequence-region directive found in "{gff_file}". '
                 f'Expected something like:\n'
                 f'"##sequence-region N1063.1 1 4417850"')
    if not genes:
        sys.exit(f'No CDS entries found in "{gff_file}".')

    return genome_len, genes


class MissingInfoFileError(Exception):
    pass


def read_nucleotide_sequences(fasta: Path, genome_folder: Path,
                              #gene_family_info = GENEFAMILYINFO,
                              initial_genome_info = INITIALGENOMEINFO) \
    -> tuple[dict[str, SeqRecord], dict[str, SeqRecord]]:
    """
    Return a dictionary mapping the gene id to its SeqRecord.
    The file `gene_family_info` contains coordinates for the genes in the root
    genome of the phylogeny. The `fasta` file is expected to contain the sequence
    of the root genome, along with the intergenic sequences

        https://biopython.org/docs/latest/api/Bio.SeqRecord.html

    Parameters
    ----------
    fasta : Path
        the file with the sequence
    genome_folder : str
        the folder containing the `gene_family_info` file
    initial_genome_info : str
        the filename containing the indices of the genes and divisions at
        the root

    Returns
    -------
    tuple[dict[str, SeqRecord], dict[str, SeqRecord]]
        map genome id (not GFF ID) to SeqRecord for the gene
        map intergene division id to SeqRecord for the gene
    """
    sequence: SeqRecord|None = None
    for i, seq_record in enumerate(SeqIO.parse(fasta, "fasta")):
        if i:
            print(f'Warning: using only the first of several entries in "{fasta}".')
            break

        sequence = seq_record

    if not sequence:
        raise(Exception(f'No sequence in "{fasta}".'))

    gidTOseq: dict[str, SeqRecord] = {}
    #gene_family_p = Path(genome_folder, gene_family_info)
    #if not gene_family_p.exists():
    #    raise MissingInfoFileError(gene_family_p)

    #print(f'getting sequence from: {gene_family_p}')
    #    #Get info about the gene families (which may have originated genes):
    #with open(gene_family_p) as f:
    #    f.readline()
    #    for line in f:
    #        gid, _, start, end = line.strip().split("\t")
    #        gidTOseq[gid] = sequence[int(start)-1: int(end)]

    init_genome_p = genome_folder / initial_genome_info
    if not init_genome_p.exists():
        raise MissingInfoFileError(init_genome_p)

        #Get info about the intergene divisions:
    didTOseq: dict[str, SeqRecord] = {}
    with open(init_genome_p) as f:
        f.readline()
        for line in f:
            type, did, start, end = line.strip().split("\t")
            start = int(start)
            end = int(end)

            if start > len(sequence)-1 or end > len(sequence):
                sys.exit(f'You have specified features in the GFF file that '
                         f'are larger that the sequence length {len(sequence)} '
                         f'(from the fasta file):\n{line}')
            
            if type == "DIVISION":
                didTOseq[did] = sequence[start-1: end]
            elif type == "GENE_FAMILY":
                gidTOseq[did] = sequence[start-1: end]
            else:
                raise(Exception(f'Unexpected TYPE "{type}" in "{init_genome_p}"'))

    return gidTOseq, didTOseq


def read_protein_sequences(gff_file: str, genome_folder: Path,
                           gene_family_info = GENEFAMILYINFO) -> dict[str, str]:
    """
    Return a dictionary mapping the gene id to its sequence.
    The file `gene_family_info` contains a the GFF IDs necessary to find
    the correct genomes in the GFF file.

    Parameters
    ----------
    gff_file : str
        the file containing features with "translated" attributes
    genome_folder : str
        the folder containing the `gene_family_info` file
    gene_family_info : str
        the filename containing the gene GFF ID and start and end indices

    Returns
    -------
    dict[str, SeqRecord]
        map genome id (not GFF ID) to SeqRecord for the gene
    """
    gffidTOseq: dict[str, str] = {}
    try:
        for rec in GFF.parse(gff_file, target_lines=1000):
            for feature in rec.features:
                if feature.type == 'CDS':
                    if 'ID' not in feature.qualifiers:
                        raise(Exception(f'Feature missing "ID" attribute.'))
                    gffid = feature.qualifiers['ID'][0]
                    if 'translation' not in feature.qualifiers:
                        raise(Exception(f'Feature missing "tranlation" attribute.'))

                    gffidTOseq[gffid] = feature.qualifiers['translation'][0]
    except TypeError as e:
        sys.exit(f'Problem opening GFF file (remove ##FASTA lines?):\n{e}')

    gidTOseq: dict[str, str] = {}
    with open(genome_folder / gene_family_info) as f:
        f.readline()
        for line in f:
            gid, gffid, _, _ = line.strip().split("\t")
            gidTOseq[gid] = gffidTOseq[gffid]

    return gidTOseq


def get_leaves_from_file(leavesfile: Path) -> set[str]:
    leaves = set()
    with open(leavesfile) as f:
        f.readline()    #Burn the opening line
        for line in f:
            leaves.add(line.strip())

    return leaves


def affected_indices_wrap(affected_indices: list[int]) -> None|tuple[int, int]:
    """
    Do the affected indices wrap around the genome end? If so, return
    (start, end) indices, where start is the start of the wrapping interval
    and end is the end (i.e. start > end).
    """
    for i, j in pairwise(sorted(affected_indices)):
        if j != i+1:
            return (j, i)

    return None


def well_behaved_indices(affected_indices: list[int]) -> bool:
    """
    We either have a series of consecutive integers, or we have two series
    of consecutive integers, where the first series is higher than the
    second.
    """
    series = 1
    for i, j in pairwise(sorted(affected_indices)):
        if j != i+1:
            series += 1
            if series > 2:
                return False

    if series == 1:
       return True

    return affected_indices[0] > affected_indices[-1]


def organize_genomes_by_branch(dir: Path, allgenomes: bool = False) \
    -> tuple[dict[str, list[tuple[int|str, Path]]], int]:
    """
    Given a directory holding genome files, return a mapping from branch
    (node) names to a list of tuples (sortkey, genome file).
    Depending on whether we are checking the "All_genomes" folder or the
    regular "Genomes" folder, the sortkey will be the sequence number, or just
    the GENOMEsuffix (every list will have length 1 in this case so the sortkey
    is not necessary).
    """
    if allgenomes:
        nodere = re.compile(rf'(\w+)-(\d+){GENOMEsuffix}$')
    else:
        nodere = re.compile(rf'(\w+)({GENOMEsuffix})$')

    node2files: dict[str, list[tuple[int|str, Path]]] = defaultdict(list)
    numfiles = 0
    for genomefile in dir.glob(f'*{GENOMEsuffix}'):
        numfiles += 1
        if m := nodere.match(genomefile.name):
            node, sortkey = m.groups()
            if sortkey != GENOMEsuffix:
                sortkey = int(sortkey)

            node2files[node].append((sortkey, genomefile))
        else:
            raise AssertionError(f'Unexpected genome file name: {genomefile.name}')

    assert numfiles, f'No GENOME files found in {dir}!'

    return node2files, numfiles


def get_last_genome(fileprefix: str) -> list[str]:
  """
  Get the genome from the filename with with the highest sequence number for the
  given path prefix.
  """
  files = glob.glob(f'{fileprefix}*{GENOMEsuffix}')
  num2file = {}
  for file in files:
    suffix = file.replace(fileprefix, '')
    num = int(suffix.replace(GENOMEsuffix, ''))
    num2file[num] = file

  return get_genome(num2file[max(num2file.keys())])


def get_genome(tsvfile: Path) -> list[str]:
  """
  Get the genome from the given filename.
  """
  #Use pandas to read the TSV file
  df = pd.read_csv(tsvfile, sep='\t')

  assert 'GENE_FAMILY' in df.columns, f'No GENE_FAMILY column in {tsvfile}!'
  assert 'ORIENTATION' in df.columns, f'No ORIENTATION column in {tsvfile}!'
  return [f'{sign}{gene}'
          for gene, sign in zip(df['GENE_FAMILY'], df['ORIENTATION'])]

