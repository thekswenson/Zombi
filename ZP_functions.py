# Load in required packages. 
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

def mkdir_cstm(path):
    # search for and, if needed, make a new directory.
    if os.path.isdir(path):
        shutil.rmtree(path)
        os.mkdir(path)
    else:
        os.mkdir(path)

def extract_files(path, pattern):
    # get every file in the directory. 
    all_files = os.listdir(path)

    # extract just the files of interest. 
    filt_files = []
    for file in all_files:
        if re.match(pattern, file):
            filt_files.append(file)

    # Sort the files.
    filt_files = sorted(filt_files)

    return filt_files

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
                    # A crude method for preventing it from removing losses. 
                    if abs(max(dists) - dists[i]) < 1:
                        # If not equal, extend branch to make it equal. 
                        node.dist = float(D(str(node.dist)) + D(str(max(dists))) - D(str(dists[i])))
                i += 1
    return tree

def apply_ultrametric(events, trees, output, rep, type = "Gene"):
    for event, tree in zip(events, trees):
    
        # Load the tree.
        t = Tree(output + f"/zombi_output/{rep}/G/{type}_trees/" + tree, format = 1)
    
        # If the tree is a single leaf (caused by origination or transfer) skip it. 
        if len(t) < 2:
            continue
    
        # Force the tree to be ultrametric. 
        t = make_ultrametric(t)
    
        # Save the tree. 
        t.write(outfile = output + f"/zombi_output/{rep}/G/sym_{type}_trees/" + tree, format = 1)
        shutil.copyfile(output + f"/zombi_output/{rep}/G/{type}_families/" + event,
                       output + f"/zombi_output/{rep}/G/sym_{type}_events/" + event)

def write_nexus(rep, output, type = "Gene"):
    # Create a dictionary to translate Zombi events to Simphy language (999 is used as a specific ignore key).
    event_keys = {"S":0, "D":1, "L":2, "T":4, "I":999, "F":999, "P":999, "O":999}
    
    # Get every tree file for that rep. 
    sym_trees = sorted(os.listdir(output + f"/zombi_output/{rep}/G/sym_{type}_trees/"))
    sym_events = sorted(os.listdir(output + f"/zombi_output/{rep}/G/sym_{type}_events/"))
    
    # Write all of the simphy prepped trees into a single file. 
    with open(output + f'/simphy_output/rep{rep}_{type}.nex','w') as wfd:
        # Create the Nexus headers 
        wfd.write('#NEXUS\n')
        wfd.write('begin trees;\n')
        for i in range(len(sym_trees)):
            # Write the tree introduction
            wfd.write(f'\ttree tree_{i+1} = [&R] ')
    
            # Load the simphy prepped tree. 
            txt = Path(output + f"/zombi_output/{rep}/G/sym_{type}_trees/" + sym_trees[i]).read_text()
    
            # Load the events file for the correct tree.
            events=pd.read_csv(output + f"/zombi_output/{rep}/G/sym_{type}_events/" + sym_events[i],sep='\t')
            
            # Now replace the node column in events with the shortened names. 
            simple_nodes = []
            for node in events.NODES:
                n = re.sub(";","_", node)
                simple_nodes.append("_".join(n.split("_")[0:2]))
            events["NODES"] = simple_nodes
    
            # Load up the U_mults for each node. 
            u_mults = pd.read_csv(output + f"/simphy_output/species_tree/{rep}/sp_rates.csv")
    
            # Get all of the nodes from the text. 
            nodes_full = re.findall('n[0-9]+_[0-9]+|int[0-9]+_[0-9]+', txt)
            nodes_trunc = re.findall('n[0-9]+|int[0-9]+', txt)
            nodes = pd.DataFrame({'NODES': nodes_full, 'nodes': nodes_trunc})
    
            # Merge the u_mult data and the event data with the IN ORDER nodes. 
            merged_1 = pd.merge(nodes, u_mults)
            merged_2 = pd.merge(merged_1, events)

            # Remove inversion and transpositions as they do not create new nodes. 
            merged_2 = merged_2.query("EVENT != 'I' and EVENT != 'P'").reset_index(drop=True)
        
            # Create a split that gets the location at the end of each branch length.
            split = re.split('([\),])', txt)

            # For each node, add in the correct parameters after the branch length. 
            cur_n = 0 
            all_params = []
            for sp in split:
                if sp == "," or sp == ")":
                    event = event_keys[merged_2.EVENT[cur_n]]
                    paralog = str(re.findall(r'[0-9]+$', merged_2.NODES[cur_n])[0])
                    u_mult = merged_2.u_mults[cur_n]
                    if event == 999:
                        new_sp = f'[&paralog={paralog},u_mult={u_mult}]{sp}' 
                    else:
                        new_sp = f'[&paralog={paralog},kind_n={event},u_mult={u_mult}]{sp}' 
                    all_params.append(new_sp)
                    cur_n += 1
                else:
                    all_params.append(sp)
            txt = ''.join(all_params)

            # Remove internal nodes. 
            txt = re.sub(r"\)int[0-9]+_[0-9]+", ")", txt)
            txt = re.sub(r"\)n[0-9]+_[0-9]+", ")", txt)
            txt = re.sub(r"_[0-9]+", "", txt)
            
            # Write to file.
            wfd.write(txt)
            wfd.write('\n')
        wfd.write('end;')

# For each rep, multiply all values in the tree by a gene specific rate modifier.
def simphy_mult(gamma_values, locus_trees_path, rep, output):
    # Find the correct alpha. 
    alpha = gamma_values["gene_alpha"][rep - 1]

    # Get the list of all trees for that rep. 
    all_s_trees = os.listdir(output + locus_trees_path + f"{rep}/1")

    for tree in all_s_trees:
        # Generate a gene specific rate modifier.
        g_mult = round(random.gammavariate(alpha,1/alpha), 4)
        if re.match("g_.*",tree):

            # Load up the tree.
            txt = Path(output + locus_trees_path + f"{rep}/1/"+ tree).read_text()

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
            with open(output + locus_trees_path + f"{rep}/"+ tree,'w') as wfd:
                wfd.write(txt)
                
    # We delete the old folder because Simphy causes a very awkward naming scheme. 
    shutil.rmtree(output + locus_trees_path + f"{rep}/1")
