# Load in required packages. 
import argparse
import subprocess
from ete3 import Tree
import zombi.AuxiliarFunctions as af
from pathlib import Path
import pandas as pd 
import re
import random
from decimal import Decimal as D
import sys
import shutil
import os
import numpy as np

# Create the argument parser. 
parser = argparse.ArgumentParser()

# Add all of the argument parser options. 
parser.add_argument("-s_loc", "--simphy_location", type=str)
parser.add_argument("-g_loc", "--genome_location", type=str)
parser.add_argument("-params", "--parameters", type=str)
parser.add_argument("-o","--output", type=str)
parser.add_argument("-fc", "--feature_choice", type=str)

args = parser.parse_args()

s_loc = args.simphy_location
g_loc = args.genome_location
params = args.parameters
output = args.output
feature_choice = args.feature_choice

def make_ultrametric(tree):
    distance = 0
    dists = []
    # Get a list of distances to the root for each leaf. 
    for post, node in tree.iter_prepostorder():
        # rounding necessary to prevent ete3 write errors down the line. 
        node.dist = round(node.dist, 1) ###!!!### This step could cause innaccuracies in simulations with few generations.
        # Create the distance list. 
        if post:
            distance -= D(str(node.dist))
        else:
            if not node.is_leaf():
                distance += D(str(node.dist))
            else:
                dists.append(float(distance + D(str(node.dist))))
    i = 0
    # Go through the tree a second time. 
    for post, node in tree.iter_prepostorder():
        if post:
            continue
        else:
            if node.is_leaf():
                # Check to make sure that the distance from the leaf to the root is equal 
                # to the max distance in our list
                if max(dists) != dists[i]:
                    # If not equal, extend branch to make it equal. 
                    node.dist = float(D(str(node.dist)) + D(str(max(dists))) - D(str(dists[i])))
                i += 1
    return tree

#################################### Generating Species Tree with Simphy ##########################################

# Configure the output directory
if os.path.isdir(output + f"/simphy_output/"):
    shutil.rmtree(output + f"/simphy_output/")
    os.mkdir(output + f"/simphy_output/")
else:
    os.mkdir(output + f"/simphy_output/")

# Read in the parameters from file and put them into a list. 
tmp_params = af.read_parameters(params,program="simphy")
reps = [i+1 for i in range(int(str(tmp_params["-RS"])))]

# Convert to a list of dictionary elements. 
s_params = []
for k,v in zip(tmp_params.keys(), tmp_params.values()):
    s_params.extend([k,v])

# Insert the necessary extra elements. 
s_params.insert(0,s_loc)
s_params.extend(["-V","3", "-O", output + "/simphy_output/species_tree"])

# Run Simphy itself using subprocess
process_txt = subprocess.run(s_params, capture_output = True, check = True)

# Save the txt file in log.txt. 
with open(output + '/simphy_output/log.txt', 'wb') as f:
    f.write(process_txt.stdout)

########################### Extract Rate Heterogeneity Parameters From SimPhy Logs ################################

# Load the Simphy log output. 
txt = Path(output + "/simphy_output/log.txt").read_text()

# extract the values from simphy's log output. 
sp_alphas = re.findall(r'-Lineage \(species\) specific rate heterogeneity gamma shape: (.*?)\n',
           txt, flags=re.DOTALL)[1:]
gene_alphas = re.findall(r'-Gene family \(locus tree\) specific rate heterogeneity gamma shape: (.*?)\n',
           txt, flags=re.DOTALL)[1:]

# Convert to numeric
sp_alphas = [eval(i) for i in sp_alphas]
gene_alphas = [eval(i) for i in gene_alphas]

# Take the values and save them as a csv.
data = {
    "rep":reps,
    "sp_alpha": sp_alphas,
    "gene_alpha": gene_alphas
}
df = pd.DataFrame(data=data)
df.to_csv(output + "/simphy_output/gamma_values.csv", index = False)

#################################### Add internal node labels ##########################################

for rep in reps:
    with open(output + f'/simphy_output/species_tree/{rep}/sz_tree.nwk','w') as wfd:
        # Load the species tree
        txt = Path(output + f"/simphy_output/species_tree/{rep}/s_tree.trees").read_text()
        
        # Load up the correct alpha value for the lineage specific gammas.
        alpha = df["sp_alpha"][rep - 1]
        
        # Create empty lists to store the node names and lineage specific rate modifier values. 
        nodes = []
        u_mults = []
        
        # Perform two splits. One to catch the number lablled leaf nodes. 
        split = re.split('([0-9]+:)', txt)
        x = 0
        for i, sp in enumerate(split):
            if not re.search('[0-9]+:', sp) is None:
                x += 1
                split[i] = re.sub('[0-9]+:', f"n{x}:",sp)
                nodes.append(f"n{x}")
                u_mults.append(round(random.gammavariate(alpha,1/alpha), 4))
        txt = ''.join(split)
        
        # And one to catch the unlabelled internal nodes. 
        split = re.split('(\):)', txt)
        x = 0
        for i, sp in enumerate(split):
            if not re.search('\):', sp) is None:
                x += 1
                split[i] = re.sub('\):', f")int{x}:",sp)
                nodes.append(f"int{x}")
                u_mults.append(round(random.gammavariate(alpha,1/alpha), 4))
        txt = ''.join(split)
        
        # save all of the internal node names and multipliers as
        data = {
            "nodes":nodes,
            "u_mults":u_mults
        }
        df2 = pd.DataFrame(data = data)
        df2.to_csv(output + f"/simphy_output/species_tree/{rep}/sp_rates.csv", index = False)
        wfd.write(txt)

#################################### Run Zombi in Genome mode for each rep ##########################################

# Configure the output directory.
if os.path.isdir(output + f"/zombi_output/"):
    shutil.rmtree(output + f"/zombi_output/")
    os.mkdir(output + f"/zombi_output/")
else:
    os.mkdir(output + f"/zombi_output/")

# Run Zombi
for rep in reps:
    to_run = ["Zombi", "Ti", output + f'/simphy_output/species_tree/{rep}/sz_tree.nwk',
              output + f'/zombi_output/{rep}']
    subprocess.run(to_run, check = True)

    to_run = ["Zombi", "Gf", params, output + f'/zombi_output/{rep}', "--genome-root",
             g_loc, "--feature-choice", feature_choice]
    subprocess.run(to_run, check = True)

################################# Force all zombi gene trees to be ultrametric ###################################

for rep in reps:
    # Get every tree file. 
    all_trees = os.listdir(output + f"/zombi_output/{rep}/G/Gene_trees")
    all_events = os.listdir(output + f"/zombi_output/{rep}/G/Gene_families/")
    all_events = sorted(all_events)
    
    # Extract just the tree files and add on the rest of the path (I must be able to improve this)
    pruned_trees = []
    for tree in all_trees:
        if re.match(r'[0-9]*_prunedtree.nwk',tree):
            pruned_trees.append(tree)
    pruned_trees = sorted(pruned_trees)

    paired_down_events = []
    for event in all_events:
        if re.match(r'[0-9]*_events.tsv',event):
            paired_down_events.append(event)
    paired_down_events = sorted(paired_down_events)
    
    # search for and, if needed, make a new directory.
    if os.path.isdir(output + f"/zombi_output/{rep}/G/sym_trees/"):
        shutil.rmtree(output + f"/zombi_output/{rep}/G/sym_trees/")
        os.mkdir(output + f"/zombi_output/{rep}/G/sym_trees/")
    else:
        os.mkdir(output + f"/zombi_output/{rep}/G/sym_trees/")
    if os.path.isdir(output + f"/zombi_output/{rep}/G/sym_events/"):
        shutil.rmtree(output + f"/zombi_output/{rep}/G/sym_events/")
        os.mkdir(output + f"/zombi_output/{rep}/G/sym_events/")
    else:
        os.mkdir(output + f"/zombi_output/{rep}/G/sym_events/")
    
    # For each of the tree files, force it to be exactly ultrametric. 
    for event, tree in zip(paired_down_events, pruned_trees):
    
        # Load the tree.
        t = Tree(output + f"/zombi_output/{rep}/G/Gene_trees/" + tree, format = 1)
    
        # If the tree is a single leaf (caused by origination or transfer) skip it. 
        if len(t) < 2:
            continue
    
        # Force the tree to be ultrametric. 
        t = make_ultrametric(t)
    
        # Save the tree. 
        t.write(outfile = output + f"/zombi_output/{rep}/G/sym_trees/" + tree, format = 1)
        shutil.copyfile(output + f"/zombi_output/{rep}/G/Gene_families/" + event,
                       output + f"/zombi_output/{rep}/G/sym_events/" + event)

########################## For each rep, combine Zombi gene trees into one nexus file #####################################

for rep in reps:
    # Get every tree file. 
    sym_trees = sorted(os.listdir(output + f"/zombi_output/{rep}/G/sym_trees/"))
    sym_events = sorted(os.listdir(output + f"/zombi_output/{rep}/G/sym_events/"))

    # Write all of the simphy prepped trees into a single file. 
    with open(output + f'/simphy_output/test{rep}.nex','w') as wfd:
        # Create the Nexus headers 
        wfd.write('#NEXUS\n')
        wfd.write('begin trees;\n')
        for i in range(len(sym_trees)):
            # Write the tree introduction
            wfd.write(f'\ttree tree_{i+1} = [&R] ')
    
            # Load the simphy prepped tree. 
            txt = Path(output + f"/zombi_output/{rep}/G/sym_trees/" + sym_trees[i]).read_text()
    
            # Load the events file for the correct tree.
            events=pd.read_csv(output + f"/zombi_output/{rep}/G/sym_events/" + sym_events[i],sep='\t')
            
            # Now replace the node column in events with the shortened names. 
            simple_nodes = []
            for node in events.NODES:
                n = re.sub(";","_", node)
                simple_nodes.append("_".join(n.split("_")[0:2]))
            events["NODES"] = simple_nodes
    
            # Load up the U_mults for each node. 
            u_mults = pd.read_csv(output + f"/simphy_output/species_tree/{rep}/sp_rates.csv")
    
            # Get all of the nodes from the text. 
            nodes = re.findall('n[0-9]+|int[0-9]+', txt)
    
            # Create a split that gets the location at then end of each branch length.
            split = re.split('([\),])', txt)
    
            # Create an integer that will count up the node we are on.
            cur_n = 0 
    
            # iterate through the split, inserting accurate U_mults, then join the split into new text.
            u_multed = []
            for sp in split:
                if sp == "," or sp == ")":
                    u_mult = u_mults[u_mults["nodes"] == nodes[cur_n]]["u_mults"].reset_index(drop = True)
                    u_multed.append(f'[&u_mult={u_mult[0]}]'+sp)
                    cur_n += 1
                else:
                    u_multed.append(sp)
            txt = ''.join(u_multed)
            
            # Replace each internal node with the type of node it is. 
            for event, node in zip(events.EVENT, events.NODES):
                if node in txt:
                    if event == "S":
                        txt = re.sub(f"{node}:", "[&kind_n=0]:", txt)
                    elif event == "D":
                        txt = re.sub(f"{node}:", "[&kind_n=1]:", txt)
                    elif event == "L":
                        txt = re.sub(f"{node}:", "[&kind_n=2]:", txt)
                    elif event == "T":
                        txt = re.sub(f"{node}:", "[&kind_n=4]:", txt)
    
            # Replace the underscore numbers with ids. 
            ids = re.findall('_[0-9]+', txt)
            for j in ids:
                txt = re.sub(j,f'[&paralog={j[1:]}]', txt)
    
            # Write to file.
            wfd.write(txt)
            wfd.write('\n')
        wfd.write('end;')

########################## Run Simphy on the Zombi gene trees #####################################

# Filter out the simphy parameters that are not needed for the second run. 
tmp_params = {key: tmp_params[key] for key in ["-SP", "-CS"]}

# Convert to a list of dictionary elements again. 
s_params = []
for k,v in zip(tmp_params.keys(), tmp_params.values()):
    s_params.extend([k,v])

# add in the necessary extra elements. 
s_params.insert(0,s_loc)
os.mkdir(output + f"/simphy_output/locus_tree")

# For each nexus file you made, run simphy using it as a locus trees input. 
for rep in reps:
    params_rep = s_params.copy()
    params_rep.extend(["-LR", output + f'/simphy_output/test{rep}.nex'])
    params_rep.extend(["-O",output + f"/simphy_output/locus_tree/{rep}"])
    subprocess.run(params_rep, check = True)

########################## Apply the lineage specific rates #####################################

# For each rep, multiply all values in the tree by a gene specific rate modifier.
for rep in reps:
    
    # Find the correct alpha. 
    alpha = df["gene_alpha"][rep - 1]

    # Get the list of all trees for that rep. 
    all_s_trees = os.listdir(output + f"/simphy_output/locus_tree/{rep}/{1}")
    
    for tree in all_s_trees:
        # Generate a gene specific rate modifier.
        g_mult = round(random.gammavariate(alpha,1/alpha), 4)
        if re.match("g_.*",tree):

            # Load up the tree.
            txt = Path(output + f"/simphy_output/locus_tree/{rep}/{1}/"+ tree).read_text()

            # Create a split to extract all values. 
            split = re.split("([0-9]+\.[0-9]+)", txt)

            # Multiply each value in the tree by the gene specific rate modifier. 
            for i, sp in enumerate(split):
                try:
                    split[i] = str(g_mult*float(sp))
                except:
                    continue
            txt = ''.join(split)

            # Save the new trees. 
            with open(output + f"/simphy_output/locus_tree/{rep}/"+ tree,'w') as wfd:
                wfd.write(txt)
    shutil.rmtree(output + f"/simphy_output/locus_tree/{rep}/{1}")