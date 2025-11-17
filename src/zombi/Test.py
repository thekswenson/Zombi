"""
Function used for testing the code.
"""
import filecmp
import pandas as pd

from collections import Counter
from enum import StrEnum
from pathlib import Path

from .AuxiliarFunctions import get_genome, get_genome_from_pieces, organize_genomes_by_branch
from .Filenames import GENOMEsuffix, PIECESsuffix

## || ## || ## || ## || ## || ## || ## || ## || ## || ## || ## || ## || ## || ##
# Comparing files


class Filetype(StrEnum):
  GENOME = 'GENOME'
  PIECES = 'PIECES'
def crosscheckGenomes(genome_folder: Path, filetype=Filetype.GENOME):
  """
  Ensure that the genomes in the `All_genomes` folder matches those in the
  `Genomes` folder.
  
  Edges (a.k.a. lineages) in the species tree are named by their pendant node
  name, and each genome on the edge has a sequential number starting from 0. If
  `n3-5` has the largest sequence number for node `n3`, then there are 6 genomes
  on the `n3` branch, and `All_genomes/n3-5_GENOME.tsv` should match
  `Genomes/n3_GENOME.tsv`.
  """
  all_genomes_folder = genome_folder / 'All_genomes'
  genomes_folder = genome_folder / 'Genomes'

  #Get the maximum number for each node:
  node2max = Counter()
  for genome_file in all_genomes_folder.glob(f'*_{filetype}.tsv'):
    base_name = genome_file.stem.replace(f'_{filetype}', '')
    node, rep = base_name.split('-')
    if node2max[node] < int(rep):
      node2max[node] = int(rep)

  #Compare files:
  for genome_file in genomes_folder.glob(f'*_{filetype}.tsv'):
    node = genome_file.stem.replace(f'_{filetype}', '')
    maxrep = node2max[node]

    all_genome_file = all_genomes_folder / f'{node}-{maxrep}_{filetype}.tsv'

    assert all_genome_file.exists(), f'Missing file: {all_genome_file}'
    assert filecmp.cmp(genome_file, all_genome_file), (
      f'{genome_file} and {all_genome_file} differ!')

    ##Ensure that there is not an extra file laying around:
    ## (this is a false alarm in the case of a transposition that puts the
    ##  transposed genes back in the same spot)
    #if maxrep > 0:
    #  previous_file = all_genomes_folder / f'{node}-{maxrep-1}_{filetype}.tsv'
    #  assert not filecmp.cmp(genome_file, previous_file), (
    #    f'{genome_file} and {previous_file} differ!')


def comparePiecesToGenomes(dir: Path, allgenomes=False):
  """
  Ensure that the gene orders are the same in the PIECES and GENOMES files.
  """
  #Depending on whether we are checking the All_genomes folder or not, the
  #sortkey will be the rep number, or just the GENOMEsuffix.
  node2files, numfiles = organize_genomes_by_branch(dir, allgenomes)

  for filelist in node2files.values():
    for _, genomefile in sorted(filelist):
      piecesname = genomefile.name.replace(GENOMEsuffix, PIECESsuffix)
      piecesfile = dir / piecesname
      assert piecesfile.exists(), f'Missing pieces file for {genomefile}'

      porder = get_genome_from_pieces(piecesfile)
      gorder = get_genome(genomefile)

      assert porder == gorder, (f'Gene order mismatch between {piecesfile} and '
                                f'{genomefile}')
      numfiles -= 1

  assert numfiles == 0

  for piecesfile in dir.glob(f'*{PIECESsuffix}'):
    genomename = piecesfile.name.replace(PIECESsuffix, GENOMEsuffix)
    assert (dir / genomename).exists(), f'Missing genome file for {piecesfile}'



#-- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - -- - - --
