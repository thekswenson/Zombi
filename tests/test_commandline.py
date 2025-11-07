"""
Tests for the command-line interface of Zombi.
We test that some of the files are created, and that the genomes created are 
consistent when using the --all-genomes flag.
"""
import pytest

from pathlib import Path

from zombi.Filenames import COMPLETETREE, TREEEVENTS, TREELENGTHS, EXTANTTREE
from zombi.Filenames import TRANSFERRATES, EVENTRATES, EXTENSIONRATES
from zombi.Test import crosscheckGenomes, comparePiecesToGenomes, Filetype

T_PARAMS = Path('Parameters/SpeciesTreeParameters.tsv')
G_PARAMS = Path('Parameters/GenomeParameters.tsv')
#S_PARAMS = Path('Parameters/SequenceParameters.tsv')
#T_PARAMS = Path('tests/SpeciesTreeParametersSeeded.tsv') #With Seed set
#G_PARAMS = Path('tests/GenomeParametersSeeded.tsv')      #With Seed set
#S_PARAMS = Path('tests/SequenceParametersSeeded.tsv')    #With Seed set
G_PARAMS_ALL = Path('tests/GenomeParametersAllgenomes.tsv')
T_SMALL_PARAMS = Path('tests/SpeciesTreeParametersSmall.tsv')
Gm_PARAMS = Path('tests/GenomeParametersAllgenomes.tsv')

@pytest.fixture(scope='session')
def projdir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project')

@pytest.fixture
def run_T(projdir, script_runner):
  """ Test the T mode of Zombi. """
  result = script_runner.run(['zombi', 'T', '-f', T_PARAMS, projdir])
  assert result.success
  return projdir / 'T'

@pytest.fixture(scope='session')
def smallprojdir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project_small')

@pytest.fixture
def small_T(smallprojdir, script_runner):
  """ Test the T mode of Zombi. """
  completetree = smallprojdir / 'T' / COMPLETETREE

  while True:
    result = script_runner.run(['zombi', 'T', '-f', T_SMALL_PARAMS, smallprojdir])
    assert result.success
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

  result = script_runner.run(['zombi', 'G', '-f', G_PARAMS_ALL, projdir])
  assert result.success

  outdir = projdir / 'G'
  assert (outdir / 'Genomes').exists()
  crosscheckGenomes(outdir)


def test_Gf(projdir, script_runner, run_T):
  """ Test the Gf mode of Zombi. """
  assert (run_T).exists(), 'There was a problem with run_T!'

  result = script_runner.run(['zombi', 'Gf', '-f', G_PARAMS_ALL, projdir])
  assert result.success

  genomedir = projdir / 'G'
  comparePiecesToGenomes(genomedir / 'Genomes')
  comparePiecesToGenomes(genomedir / 'All_genomes', True)
  crosscheckGenomes(genomedir)
  crosscheckGenomes(genomedir, Filetype.PIECES)


@pytest.fixture
def run_RateCustomizer(projdir, script_runner):
  """ Test the RateCustomizer mode of Zombi. """
  result = script_runner.run(['zombiRateCustomizer', 'G', G_PARAMS, projdir])
  assert result.success

  customrates = projdir / 'CustomRates'
  assert (customrates / TRANSFERRATES).exists()
  assert (customrates / EVENTRATES).exists()
  assert (customrates / EXTENSIONRATES).exists()
  return True


def test_Gu(projdir, script_runner, run_T, run_RateCustomizer):
  """ Test the Gu mode of Zombi. """
  assert run_RateCustomizer, 'There was a problem with run_RateCustomizer!'
  assert (run_T).exists(), 'There was a problem with run_T!'

  result = script_runner.run(['zombi', 'Gu', '-f', G_PARAMS_ALL, projdir])
  assert result.success

  outdir = projdir / 'G'
  assert (outdir / 'Genomes').exists()
  crosscheckGenomes(outdir)


def test_Gm(smallprojdir, script_runner, small_T):
  """ Test the Gm mode of Zombi. """
  assert (small_T).exists(), 'There was a problem with run_T!'

  result = script_runner.run(['zombi', 'Gm', '-f', Gm_PARAMS, smallprojdir])
  assert result.success

  outdir = smallprojdir / 'G'
  assert (outdir / 'Genomes').exists()
  #crosscheckGenomes(outdir)



#_______________________________________________________________________________
# Functions
