import os
import re
import argparse
import subprocess 
import shutil
from pathlib import Path

# Create the argument parser. 
parser = argparse.ArgumentParser()

# Add all of the argument parser options. 
parser.add_argument("-g_loc", "--genome_location", type=str)
parser.add_argument("-params", "--parameters", type=str)
parser.add_argument("-o","--output", type=str)
parser.add_argument("-r", "--reps", type=str)
parser.add_argument("-p","--parallel", type=str)

args = parser.parse_args()

seq_root = args.genome_location
seq_params = args.parameters
output = args.output
reps = args.reps
parallel = args.parallel
reps = [i+1 for i in range(int(reps))]

# Run Zombi in S mode for all reps. 
for rep in reps:
    print(f"##### Running Sequence Simulation for rep {rep} #####")
    to_run = ["Zombi", "Sf", seq_params, output + f'/zombi_output/{rep}',
              "--sequence-root", seq_root, "--parallel", parallel]
    subprocess.run(to_run, capture_output = True, check = True)

print("# Done")