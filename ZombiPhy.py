# Load in required packages. 
import argparse
import subprocess
from ete3 import Tree
import zombi.AuxiliarFunctions as af
import ZP_functions as zp
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
gamma_values = pd.DataFrame(data=data)
gamma_values.to_csv(output + "/simphy_output/gamma_values.csv", index = False)

#################################### Add internal node labels ##########################################

for rep in reps:
    with open(output + f'/simphy_output/species_tree/{rep}/sz_tree.nwk','w') as wfd:
        # Load the species tree
        txt = Path(output + f"/simphy_output/species_tree/{rep}/s_tree.trees").read_text()
        
        # Load up the correct alpha value for the lineage specific gammas.
        alpha = gamma_values["sp_alpha"][rep - 1]
        
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
        df = pd.DataFrame(data = data)
        df.to_csv(output + f"/simphy_output/species_tree/{rep}/sp_rates.csv", index = False)
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

    # Extract lists of all of the relevant files. 
    all_gene_trees = zp.extract_files(output + f"/zombi_output/{rep}/G/Gene_trees", r'[0-9]*_prunedtree.nwk')
    all_division_trees = zp.extract_files(output + f"/zombi_output/{rep}/G/Division_trees", r'[0-9]*_prunedtree.nwk')
    all_gene_events = zp.extract_files(output + f"/zombi_output/{rep}/G/Gene_families", r'[0-9]*_events.tsv')
    all_division_events = zp.extract_files(output + f"/zombi_output/{rep}/G/Division_families", r'[0-9]*_events.tsv')
    
    # Make new directories to hold the edited trees and events. 
    zp.mkdir_cstm(output + f"/zombi_output/{rep}/G/sym_Gene_trees/")
    zp.mkdir_cstm(output + f"/zombi_output/{rep}/G/sym_Gene_events/")
    zp.mkdir_cstm(output + f"/zombi_output/{rep}/G/sym_Division_trees/")
    zp.mkdir_cstm(output + f"/zombi_output/{rep}/G/sym_Division_events/")

    # Force the trees to be ultrametric.
    zp.apply_ultrametric(all_gene_events, all_gene_trees, output, rep, "Gene")
    zp.apply_ultrametric(all_division_events, all_division_trees, output, rep, "Division")

########################## For each rep, combine Zombi gene trees into one nexus file #####################################

for rep in reps:
    zp.write_nexus(rep, output, "Gene")
    zp.write_nexus(rep, output, "Division")

########################## Run Simphy on the Zombi gene trees #####################################

# Filter out the simphy parameters that are not needed for the second run. 
tmp_params = {key: tmp_params[key] for key in ["-SP", "-CS"]}

# Convert to a list of dictionary elements again. 
s_params = []
for k,v in zip(tmp_params.keys(), tmp_params.values()):
    s_params.extend([k,v])

# add in the necessary extra elements. 
s_params.insert(0,s_loc)

# Configure the output directory
zp.mkdir_cstm(output + f"/simphy_output/locus_tree")
zp.mkdir_cstm(output + f"/simphy_output/locus_tree/Gene")
zp.mkdir_cstm(output + f"/simphy_output/locus_tree/Division")

# For each nexus file you made, run simphy using it as a locus trees input. 
for rep in reps:
    params_rep = s_params.copy()
    params_rep.extend(["-LR", output + f'/simphy_output/rep{rep}_Gene.nex'])
    params_rep.extend(["-O",output + f"/simphy_output/locus_tree/Gene/{rep}"])
    subprocess.run(params_rep, check = True)
    params_rep = s_params.copy()
    params_rep.extend(["-LR", output + f'/simphy_output/rep{rep}_Division.nex'])
    params_rep.extend(["-O",output + f"/simphy_output/locus_tree/Division/{rep}"])
    subprocess.run(params_rep, check = True)

########################## Apply the lineage specific rates #####################################
    
for rep in reps:
    zp.simphy_mult(gamma_values, "/simphy_output/locus_tree/Gene/", rep, output)
    zp.simphy_mult(gamma_values, "/simphy_output/locus_tree/Division/", rep, output)