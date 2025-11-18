#!/usr/bin/env python3
"""
Function to analyze and compare the positional orthologs computed by the
various methods.
"""
import random
import sys
import argparse
import json
import re

from collections import defaultdict
from itertools import combinations
from pathlib import Path

import pandas as pd
from sklearn.metrics import rand_score


def mapMethodToFilename(files: list[str]) -> dict[str, str]:
  """
  Map the method name to the filename.
  """
  method2file = {}
  for f in files:
    method = f.split('/')[1]
    method2file[method] = f

  return method2file


def makeCopynumToRandFiles(duplist: list[str], randlist: list[str],
                           outfile: str) -> None:
  """
  Given a list of duplication files and a list of random files, create a
  dictionary mapping the copy number to the family ID and the Rand index.
  """
  copynum2rands = mapCopynumToRand(duplist, randlist)

  #Write the output file
  with open(outfile, 'w') as f:
    f.write('COPY_NUMBER\tFAMILY\tRAND_INDEX\n')
    for copynumber in sorted(copynum2rands.keys(), reverse=True):
      for family, rand in copynum2rands[copynumber]:
        f.write(f'{copynumber}\t{family}\t{rand}\n')


def mapCopynumToRand(duplist: list[str], randlist: list[str]) \
  -> dict[int, list[tuple[str, float]]]:
  """
  From a duplication TSV file with lines of the form:
    FAMILY_ID   TOTAL_DUPS  DUPS    TANDEM_DUPS
  and a list of TSV files of the form:
    FAMILY_ID RAND_INDEX
  """
  method2dupfile = mapMethodToFilename(duplist)
  method2randfile = mapMethodToFilename(randlist)

  cn2rand = defaultdict(list)
  for method, dupfile in method2dupfile.items():
    randfile = method2randfile[method]

    fam2cn = {}                     #map family to copy number
    with open(dupfile) as f:
      f.readline()                  #header
      for line in f:
        fam, total, dups, tdups  = line.strip().split('\t')
        fam2cn[fam] = int(total)+1  #transform number of dups to copy number

    with open(randfile) as f:
      for line in f:
        fam, rand = line.strip().split('\t')

        if fam not in fam2cn:
          raise ValueError(f'Family {fam} not found in {dupfile}!')

        if rand != 'None':
          cn2rand[fam2cn[fam]].append((fam, float(rand)))

  return cn2rand


def randIndicesBetweenFiles(file1: str, file2: str, outf: str,
                            log: Path|None = None) -> None:
  """
  For all the partitions given the JSON files, compute the Rand index between
  the partitions of the same name.
  """
  #Read the JSON files:
  with open(file1) as f:
    partitions1 = json.load(f)['data'] #map partition name to partition
  with open(file2) as f:
    partitions2 = json.load(f)['data'] #map partition name to partition

  #Compute the Rand index between the partitions:
  name2rand = {}
  warnings = []
  for name in partitions1.keys():
    if name not in partitions2:
      #raise ValueError(f'Partition {name} not found in {file2}!')
      if partitions1[name]:
        warnings.append(f'WARNING: Partition {name} not found in {file2}!')

      name2rand[name] = None
    else:
      name2rand[name] = randIndex(toSet(partitions1[name]),
                                  toSet(partitions2[name]),
                                  warnings)

  #Write the Rand indices to the TSV output file:
  with open(outf, 'w') as f:
    for name, rand in name2rand.items():
      f.write(f'{name}\t{rand}\n')

  #Write the warnings to the log file:
  if log:
    with open(log, 'w') as f:
      for w in warnings:
        f.write(f'{w}\n')


def toSet(partition: list[list[list[str]]]) -> list[set[tuple[str, str]]]:
  """ Convert the parts to sets. """
  return [{(e[0], e[1]) for e in s} for s in partition]


def toSets(partitions: dict[str, list[list[list[str]]]]) \
  -> dict[str, list[set[tuple[str, str]]]]:
  """ Convert the parts to sets. """
  return {name: toSet(part) for name, part in partitions.items()}


def randIndex(p1: list[set[tuple[str, str]]],
              p2: list[set[tuple[str, str]]],
              warnings: list[str]|None=None) -> float:
  """
  Compute the Rand index between two partitions.

  The Rand index is the number of pairs of elements that are assigned to the
  same set in both partitions, plus the number of pairs of elements that are
  assigned to different sets in both partitions, divided by the total number
  of pairs of elements.
  """
  #Create single element sets if there are elements missing from a partition:
  if p1:
    p1alphabet = set.union(*p1)
  else:
    p1alphabet = set()
  if p2:
    p2alphabet = set.union(*p2)
  else:
    p2alphabet = set()
  universe = p1alphabet | p2alphabet

  if len(universe) <= 1:
    return 1.0

  p1missing = universe - p1alphabet
  if p1missing and warnings != None:
    warnings.append(f'WARNING: Partition 1 missing elements: {p1missing}')

  p2missing = universe - p2alphabet
  if p2missing and warnings != None:
    warnings.append(f'WARNING: Partition 2 missing elements: {p2missing}')

  for e in p1missing:
    p1.append(set([e]))
  for e in p2missing:
    p2.append(set([e]))
    
    
  #Set up vectors for scikit-learn to compute the Rand index:
  v1 = [0] * len(universe)
  v2 = [0] * len(universe)
  e2i = {e: i for i, e in enumerate(universe)}

  for j, s in enumerate(p1):
    for e in s:
      v1[e2i[e]] = j
  for j, s in enumerate(p2):
    for e in s:
      v2[e2i[e]] = j
    
  #Map element to set:
  e2s1 = {e: s for s in p1 for e in s}
  e2s2 = {e: s for s in p2 for e in s}

  shared = 0      #Pair in same set in both partitions
  notshared = 0   #Pair in different sets in both partitions
  bad = 0
  #Classify the pairs of elements as shared both or not shared in both:
  for e1, e2 in combinations(universe, 2):
    if e2s1[e1] == e2s1[e2]:
      if e2s2[e1] == e2s2[e2]:
        shared += 1
      else:
        bad += 1

    else:
      if e2s2[e1] == e2s2[e2]:
        bad += 1
      else:
        notshared += 1

  assert shared + notshared + bad == numPairs(len(universe))
  assert rand_score(v1, v2) == (shared + notshared) / (shared + notshared + bad)
  return (shared + notshared) / (shared + notshared + bad)

    
def numPairs(n) -> int:
  """
  Compute the number of pairs of genes in a set of n genes.
  The number of pairs is given by n * (n - 1) / 2.
  """
  return (n * (n - 1)) / 2


def createRandomPartition(inpartitions: Path, outpartitions: Path) -> None:
  """
  Given a JSON file with the true partition, make a random partition containing
  the same number of parts.  A valid partition has at most one element from
  each genome in each part.  The random partition will have the same number of
  parts as the true partition.

  Each element is a tuple of the form (genome, gene).
  """
  #Read the JSON
  with open(inpartitions) as f:
    partitions = json.load(f)['data'] #map partition name to partition
  
  newpartitions = {}
  for name, partition in partitions.items():
    if not partition:
      print(f'WARNING: Partition {name} is empty!')
      newpartitions[name] = []
      continue

    #Organize elements according to their genome:
    genome2elements = defaultdict(list)
    for part in partition:
      for e in part:
        genome2elements[e[0]].append(e)
    
    #Create the buckets:
    buckets = [[] for _ in range(len(partition))]
    #Add each element to a random bucket:
    for elements in genome2elements.values():
      #Choose len(elements) random buckets:
      random.shuffle(elements)
      for bucket, element in zip(random.sample(buckets, len(elements)), elements):
        bucket.append(element)

    newpartitions[name] = buckets

  #Write the new partitions to the output file:
  with open(outpartitions, 'w') as f:
    json.dump({'comment': 'Random partitions made by createRandomParitition()',
               'data': newpartitions}, f, indent=2)


cpnumfilere = re.compile(r'.+copynumberTOfamilyTOrand_(\S+)\((\S+)\).tsv')
def combineTablesAddMethod(tablefiles: list[str], outfile: str) -> None:
  """
  Given a list of files with the following columns:
    COPY_NUMBER FAMILY_ID RAND_INDEX
  Combine the files into a single file with the following columns:
    COPY_NUMBER FAMILY_ID RAND_INDEX METHOD
  """
  newdf = pd.DataFrame()
  #Aggregate the Rand index values from the different methods, adding a
  #column for the method:
  for f in tablefiles:
    #Get the method name from the filename:
    m = cpnumfilere.match(f)
    if m is None:
      raise ValueError(f'Could not parse {f} for method name')

    method = m.group(1)
    params = m.group(2)
    #Read in the file and add a column for the method:
    df = pd.read_csv(f, sep='\t')
    df['METHOD'] = f'{method}({params})'
    #Add the data to the new dataframe:
    newdf = pd.concat([newdf, df], ignore_index=True)

  #Write the new dataframe to a file:
  newdf.to_csv(outfile, sep='\t', index=False)


def getMatchingSimulationParameters(filename: str,
                                    paramstrings: list[str]) -> str:
  """
  Given a filename, return the first parameter string that occurs in it.
  """
  for pstring in paramstrings:
    if pstring in filename:
      return pstring

  raise ValueError(f'Could not find parameter string in {filename}.')


#PARAMETERS:    __    __    __    __    __    __    __    __    __    __    __ 
#__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \_

def main():
    desc = 'Analyze positional ortholog files.'
    parser = argparse.ArgumentParser(description=desc)
    subparsers = parser.add_subparsers(title='subcommands')

    parent = argparse.ArgumentParser(add_help=False)

    p_ri = subparsers.add_parser('randindex', aliases=['ri'], parents=[parent],
                                 help='Compute the Rand index.')
  
    p_rp = subparsers.add_parser('randomizepartition', aliases=['rp'], parents=[parent],
                                 help='Randomize the given (simulated) partition.')
 
    p_ri.add_argument('IN_FILE', type=Path, nargs=2,
                      help='input JSON')
    p_ri.add_argument('OUT_TSV', type=Path,
                      help='TSV output filename')

    p_rp.add_argument('IN_FILE', type=Path,
                      help='input JSON')
    p_rp.add_argument('OUT_FILE', type=Path,
                      help='output JSON')
 
 
    p_ri.set_defaults(func=rand_indices)
    p_rp.set_defaults(func=randomize_partition)

    args = parser.parse_args()
    
    if 'func' not in vars(args):
      parser.print_help()
      sys.exit()

    args.func(args) 


#MAIN:    __    __    __    __    __    __    __    __    __    __    __    __
#__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \_
 

def rand_indices(args):
  randIndicesBetweenFiles(args.IN_FILE[0], args.IN_FILE[1], args.OUT_TSV)

def randomize_partition(args):
  createRandomPartition(args.IN_FILE, args.OUT_FILE)
    
if __name__ == '__main__':
  main()
