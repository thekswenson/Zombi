"""
Tests for the command-line interface of Zombi.
We test that some of the files are created, and that the genomes created are 
consistent when using the --all-genomes flag.
"""
import filecmp
import pytest

from pathlib import Path
from collections import Counter
from enum import StrEnum

from zombi.Filenames import COMPLETETREE, TREEEVENTS, TREELENGTHS, EXTANTTREE
from zombi.Filenames import TRANSFERRATES, EVENTRATES, EXTENSIONRATES

T_PARAMS = Path('Parameters/SpeciesTreeParameters.tsv')
G_PARAMS = Path('Parameters/GenomeParameters.tsv')
#S_PARAMS = Path('Parameters/SequenceParameters.tsv')
#T_PARAMS = Path('tests/SpeciesTreeParameters.tsv') #With Seed set
#G_PARAMS = Path('tests/GenomeParameters.tsv')      #With Seed set
#S_PARAMS = Path('tests/SequenceParameters.tsv')    #With Seed set
T_SMALL_PARAMS = Path('tests/SpeciesTreeParameters_small.tsv')
Gm_PARAMS = Path('tests/GenomeParameters.tsv')

@pytest.fixture(scope='session')
def projdir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project')

@pytest.fixture
def run_T(projdir, script_runner):
  """ Test the T mode of Zombi. """
  script_runner.run(['zombi', 'T', T_PARAMS, projdir])
  return projdir / 'T'

@pytest.fixture(scope='session')
def smallprojdir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project_small')

@pytest.fixture
def small_T(smallprojdir, script_runner):
  """ Test the T mode of Zombi. """
  completetree = smallprojdir / 'T' / COMPLETETREE

  while True:
    script_runner.run(['zombi', 'T', '-f', T_SMALL_PARAMS, smallprojdir])
    assert completetree.exists()

    if completetree.stat().st_size > 20:
      return smallprojdir / 'T'


def test_T(run_T):
  assert (run_T / TREEEVENTS).exists()
  assert (run_T / COMPLETETREE).exists()
  assert (run_T / EXTANTTREE).exists()
  assert (run_T / TREELENGTHS).exists()


def test_G(projdir, script_runner, run_T):
  """ Test the G mode of Zombi. """
  assert (run_T).exists(), 'There was a problem with run_T!'

  script_runner.run(['zombi', 'G', '-a', G_PARAMS, projdir])

  outdir = projdir / 'G'
  assert (outdir / 'Genomes').exists()
  crosscheckGenomes(outdir)


def test_Gf(projdir, script_runner, run_T):
  """ Test the Gf mode of Zombi. """
  assert (run_T).exists(), 'There was a problem with run_T!'

  script_runner.run(['zombi', 'Gf', '-fa', G_PARAMS, projdir])

  outdir = projdir / 'G'
  assert (outdir / 'All_genomes').exists()
  crosscheckGenomes(outdir)
  crosscheckGenomes(outdir, Filetype.PIECES)


@pytest.fixture
def run_RateCustomizer(projdir, script_runner):
  """ Test the RateCustomizer mode of Zombi. """
  script_runner.run(['zombiRateCustomizer', 'G', G_PARAMS, projdir])

  customrates = projdir / 'CustomRates'
  assert (customrates / TRANSFERRATES).exists()
  assert (customrates / EVENTRATES).exists()
  assert (customrates / EXTENSIONRATES).exists()
  return True


def test_Gu(projdir, script_runner, run_T, run_RateCustomizer):
  """ Test the Gu mode of Zombi. """
  assert run_RateCustomizer, 'There was a problem with run_RateCustomizer!'
  assert (run_T).exists(), 'There was a problem with run_T!'

  script_runner.run(['zombi', 'Gu', '-fa', G_PARAMS, projdir])

  outdir = projdir / 'G'
  assert (outdir / 'Genomes').exists()
  crosscheckGenomes(outdir)


def test_Gm(smallprojdir, script_runner, small_T):
  """ Test the Gm mode of Zombi. """
  assert (small_T).exists(), 'There was a problem with run_T!'

  script_runner.run(['zombi', 'Gm', '-fa', Gm_PARAMS, smallprojdir])

  outdir = smallprojdir / 'G'
  assert (outdir / 'Genomes').exists()
  #crosscheckGenomes(outdir)



#_______________________________________________________________________________
# Functions

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

    assert all_genome_file.exists()
    assert filecmp.cmp(genome_file, all_genome_file)

    if maxrep > 0:
      previous_file = all_genomes_folder / f'{node}-{maxrep-1}_{filetype}.tsv'
      assert not filecmp.cmp(genome_file, previous_file)
