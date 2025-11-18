#!/usr/bin/env python3
"""
Make plots.
"""
import argparse
import re

from pathlib import Path
import pandas as pd

import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


def plotRandMeansByDuplication(table: Path, outfile: Path):
  """
  Given a TSV file with the following columns:
    COPY_NUMBER FAMILY RAND_INDEX METHOD
  plot the mean values of the rand index for each copy number, for each method.
  """
  #Read the dataframe:
  df = pd.read_csv(table, sep='\t')
  
  # Plot the means values for each copy number, an independent line for each method:
  sns.set(style="whitegrid")
  fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

  # Assign a different marker for each method
  methods = df['METHOD'].unique()
  markers = ['o', 's', 'D', '^', 'v', 'P', '*', 'X', '<', '>']  # extend if needed
  marker_dict = {method: markers[i % len(markers)] for i, method in enumerate(methods)}

  # Top plot: Rand Index means by copy number for each method
  for method in methods:
    method_df = df[df['METHOD'] == method]
    sns.lineplot(
      x='COPY_NUMBER',
      y='RAND_INDEX',
      data=method_df,
      label=method,
      marker=marker_dict[method],
      lw=2,
      errorbar='sd',
      #palette='Set1',
      markersize=8,
      ax=axes[0]
    )

  axes[0].set_ylabel('Rand Index')
  axes[0].set_ylim(0.0, 1.1)
  axes[0].legend(title='Method', bbox_to_anchor=(1, 1))
  #axes[0].set_title('Rand Index by Copy Number')

  # Bottom plot: Number of families per copy number for the first method
  first_method = methods[0]
  first_method_df = df[df['METHOD'] == first_method]
  family_counts = first_method_df.groupby('COPY_NUMBER')['FAMILY'].nunique().reset_index()
  sns.barplot(
    x='COPY_NUMBER',
    y='FAMILY',
    data=family_counts,
    color='gray',
    ax=axes[1]
  )
  axes[1].set_ylabel('Number of Families')
  axes[1].set_xlabel('Copy Number')

  max_copy = int(df['COPY_NUMBER'].max())
  axes[1].set_xlim(-0.5, max_copy + 0.5)
  axes[1].set_xticks(range(0, max_copy + 1))

  plt.tight_layout()
  plt.savefig(outfile)
  plt.close()


#Filename looks like: copynumberTOfamilyTOrand_{method}({spec})-{seqparams}.tsv'
cnfilere = re.compile(r'copynumberTOfamilyTOrand_(\S+)\((\S+).tsv')
def boxplotRandByDuplication(copynumberfile: Path, outfile: Path):
  """
  Given a TSV file with the following columns:
    COPY_NUMBER FAMILY_ID RAND_INDEX
  output a boxplot with the copy number on the x-axis and the rand index on the
  y-axis.
  """
  df = pd.read_csv(copynumberfile, sep='\t')

  #Build the title from the filename:
  if m := cnfilere.match(Path(copynumberfile).name):
    method = m.group(1)
    parameters = m.group(2)
    title = f'{method} ({parameters}'
  else:
    print(f'WARNING: could not parse filename {copynumberfile.name}.')
    title = 'Rand Index by Copy Number'

  #Make the boxplot using seaborn:
  sns.set(style="whitegrid")
  plt.figure(figsize=(10, 6))
  
  #Make boxplot such that y-axis shows a range of 0.5 to 1.0, unless there is
  #mean value (per copy number) that is less than 0.5, then show the range of
  #0.0 to 1.0:
  minmean = df.groupby('COPY_NUMBER')['RAND_INDEX'].mean().min()
  if minmean < 0.5:
    plt.ylim(0.0, 1.0)
  else:
    plt.ylim(0.5, 1.0)
  
  #Make boxplot:
  sns.boxplot(x='COPY_NUMBER', y='RAND_INDEX', data=df,
              showmeans=True,
              meanprops={"marker": "o", "markerfacecolor": "white",
                         "markeredgecolor": "black", "markersize": 10,
                         "linestyle": "None"},
              boxprops={"facecolor": "white", "edgecolor": "black"},
              medianprops={"color": "black", 'linewidth': 2},
              whiskerprops={"color": "black", 'linewidth': 2},
              capprops={"color": "black", 'linewidth': 2},
              showfliers=False, notch=True)

  plt.title(title)
  plt.ylabel('Rand Index')
  plt.xlabel('Copy Number')
  #plt.xticks(rotation=70)

  plt.tight_layout()

  plt.savefig(outfile)
  plt.close()


def makeRandBoxplots(files: list[str], outfile: str, labels: list[str],
                     title='') -> None:
  """
  Given a list of files mapping family id to rand index, make a boxplot.
  """
  assert len(files) == len(labels)
  #Ensure that the output directory exists:
  outpath = Path(outfile)
  outpath.parent.mkdir(parents=True, exist_ok=True)
  
  #Each file is a TSV file with two columns: family id and rand index.
  #Read each file as a vector of rand floats:
  file2randvec = {}
  for f in files:
    with open(f) as fh:
      randvec = []
      for line in fh:
        fields = line.strip().split()
        if fields[1] != 'None':
          randvec.append(float(fields[1]))

      file2randvec[f] = randvec
  
  
  #Use seaborn to make the boxplot:
  sns.set(style="whitegrid")
  plt.figure(figsize=(10, 6))

  #Make the boxplot with wedges (confidence intervals for the median):
  sns.boxplot(data=file2randvec, showmeans=True,
              meanprops={"marker": "o", "markerfacecolor": "white",
                         "markeredgecolor": "black", "markersize": 10,
                         "linestyle": "None"},
              boxprops={"facecolor": "white", "edgecolor": "black"},
              medianprops={"color": "black", 'linewidth': 2},
              whiskerprops={"color": "black", 'linewidth': 2},
              capprops={"color": "black", 'linewidth': 2},
              showfliers=False, notch=True)

  #sns.boxplot(data=file2randvec, palette="Set3")

  plt.title(title)
  plt.ylabel('Rand Index')

  #Add teh x-axis labels:
  plt.xticks(range(len(labels)), labels, rotation=70)
  #plt.xlabel('Methods')

  plt.tight_layout()

  plt.savefig(outpath)
  plt.close()

  
  


#PARAMETERS:    __    __    __    __    __    __    __    __    __    __    __ 
#__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \_

def main():
    pass
    desc = 'Make plots.'
    parser = argparse.ArgumentParser(description=desc)
    subparsers = parser.add_subparsers(title='subcommands')

    parent = argparse.ArgumentParser(add_help=False)

    p_rm = subparsers.add_parser('RandMeans', aliases=['rm'], parents=[parent],
                                 help='Compute the Rand index.')
 
    p_rm.add_argument('IN_FILE', type=Path,
                      help='input TSV')
    p_rm.add_argument('OUT_FILE', type=Path,
                      help='PDF output file')
 
    p_rm.set_defaults(func=plot_Rand_means)

    #args = parser.parse_args()
    
    #if 'func' not in vars(args):
    #  parser.print_help()
    #  sys.exit()

    #args.func(args) 


#MAIN:    __    __    __    __    __    __    __    __    __    __    __    __
#__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \_
 
def plot_Rand_means(args):
  plotRandMeansByDuplication(args.IN_FILE, args.OUT_FILE)

if __name__ == '__main__':
  main()
