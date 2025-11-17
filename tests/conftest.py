"""
The fixtures to be shared across multiple test modules.
"""
import subprocess
import pytest

from typing import NamedTuple
from pathlib import Path

from zombi.Filenames import COMPLETETREE, EVENTRATES, EXTENSIONRATES
from zombi.Filenames import TRANSFERRATES


# . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ .
# Files:
T_PARAMS = Path('Parameters/SpeciesTreeParameters.tsv')
T_SMALL_PARAMS = Path('tests/SpeciesTreeParametersSmall.tsv')
G_PARAMS = Path('Parameters/GenomeParameters.tsv')
G_PARAMS_ALL = Path('tests/GenomeParametersAllgenomes.tsv')
G_PARAMS_SEEDED = Path('tests/GenomeParametersSeeded.tsv')
S_PARAMS = Path('Parameters/SequenceParameters.tsv')


# . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ . _ .
# Fixtures:

@pytest.fixture(scope='session')
def projdir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project')

@pytest.fixture(scope='session')
def smallprojdir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project_small')


#     .     .     .     .     .     .     .     .     .     .     .     .
# Trees:

def run_T_factory(pdir, runner):
  """ Run the T mode of zombi with one of the given runners. """
  completetree = pdir / 'T' / COMPLETETREE

  result = runner.run(['zombi', 'T', T_PARAMS, pdir])
  if hasattr(result, 'success'):
    assert result.success, f'Error running T:\n{result.stderr}\n{result.stdout}'
  else:
    assert result.returncode == 0, (
      f'Error running T:\n{result.stderr}\n{result.stdout}')

  assert completetree.exists()
  return pdir / 'T'


@pytest.fixture(scope='session')
def run_T(projdir):
  """ Run the T mode of Zombi. Module-scoped to avoid rerunning. """
  return run_T_factory(projdir, subprocess)

@pytest.fixture
def rerun_T(projdir, script_runner):
  """ For a rerun of the T mode of Zombi. """
  return run_T_factory(projdir, script_runner)


def run_small_T_factory(pdir, runner):
  """ Run the T mode of zombi with one of the given runners. """
  completetree = pdir / 'T' / COMPLETETREE

  numtries = 100
  while True:
    result = runner.run(['zombi', 'T', '-f', T_SMALL_PARAMS, pdir])
    if hasattr(result, 'success'):
      assert result.success, f'Error running T:\n{result.stderr}\n{result.stdout}'
    else:
      assert result.returncode == 0, (
        f'Error running T:\n{result.stderr}\n{result.stdout}')

    assert completetree.exists()

    if completetree.stat().st_size > 20:
      return pdir / 'T'

    print(f'INFO: Tree too small ({completetree.stat().st_size} bytes), rerunning...')
    assert numtries != 0, 'Could not generate a small tree after many attempts!'
    numtries -= 1


@pytest.fixture(scope='session')
def run_small_T(smallprojdir):
  """ Test the T mode of Zombi. """
  return run_small_T_factory(smallprojdir, subprocess)

@pytest.fixture
def rerun_small_T(smallprojdir, script_runner):
  """ Test the T mode of Zombi. """
  return run_small_T_factory(smallprojdir, script_runner)


#     .     .     .     .     .     .     .     .     .     .     .     .
# Genomes:

class RunDirs(NamedTuple):
  T: Path
  G: Path

def G_runner(pdir, mode='G', force=True, allgenomes=False, seeded=False) -> Path:
  """
  Run the G mode of Zombi.

  Parameters
  ----------
  pdir : Path
      The project directory.
  mode : str, optional
      The mode to run (one of G, Gf, Gm, or Gu).
  force : bool, optional
      Whether to force overwrite existing files (True) or not (False).
  allgenomes : bool, optional
      Whether to generate all genomes (True) or just those at the nodes of
      the species tree.
  seeded : bool, optional
      Whether to seed the random number generators (True) or not (False).
  """
  params = G_PARAMS
  if allgenomes:
    params = G_PARAMS_ALL
  elif seeded:
    params = G_PARAMS_SEEDED
  fflag = ''
  if force:
    fflag = '-f'
  result = subprocess.run(['zombi', mode, fflag, params, pdir])
  assert result.returncode == 0, (
    f'Error running G:\n{result.stderr}\n{result.stdout}')

  assert (pdir / 'G').exists()
  return pdir / 'G'


@pytest.fixture(scope='session')
def run_G(projdir, run_T) -> RunDirs:
  """ Run the G mode of Zombi. Module-scoped to avoid rerunning. """
  return RunDirs(run_T, G_runner(projdir))

@pytest.fixture(scope='session')
def run_G_all(projdir, run_T) -> RunDirs:
  """ Run the G mode of Zombi. Module-scoped to avoid rerunning. """
  return RunDirs(run_T, G_runner(projdir, allgenomes=True))

@pytest.fixture(scope='session')
def run_G_seeded(projdir, run_T) -> RunDirs:
  """ Run the G mode of Zombi. Module-scoped to avoid rerunning. """
  return RunDirs(run_T, G_runner(projdir, seeded=True))

@pytest.fixture
def rerun_G(projdir, run_T) -> RunDirs:
  """ Run the G mode of Zombi. """
  return RunDirs(run_T, G_runner(projdir))

@pytest.fixture
def rerun_G_all(projdir, run_T) -> RunDirs:
  """ Run the G mode of Zombi. """
  return RunDirs(run_T, G_runner(projdir, allgenomes=True))

@pytest.fixture
def rerun_G_seeded(projdir, run_T) -> RunDirs:
  """ Run the G mode of Zombi. """
  return RunDirs(run_T, G_runner(projdir, seeded=True))


@pytest.fixture(scope='session')
def run_Gf(projdir, run_T) -> RunDirs:
  """ Run the Gf mode of Zombi. Module-scoped to avoid rerunning. """
  return RunDirs(run_T, G_runner(projdir, mode='Gf'))

@pytest.fixture(scope='session')
def run_Gf_all(projdir, run_T) -> RunDirs:
  """ Run the Gf mode of Zombi. Module-scoped to avoid rerunning. """
  return RunDirs(run_T, G_runner(projdir, mode='Gf', allgenomes=True))

@pytest.fixture
def rerun_Gf(projdir, run_T) -> RunDirs:
  """ Run the Gf mode of Zombi. """
  return RunDirs(run_T, G_runner(projdir, mode='Gf'))

@pytest.fixture
def rerun_Gf_all(projdir, run_T) -> RunDirs:
  """ Run the Gf mode of Zombi. """
  return RunDirs(run_T, G_runner(projdir, mode='Gf', allgenomes=True))


@pytest.fixture(scope='session')
def run_Gm(smallprojdir, run_small_T) -> RunDirs:
  """ Run the Gm mode of Zombi. Module-scoped to avoid rerunning. """
  return RunDirs(run_small_T, G_runner(smallprojdir, mode='Gm', allgenomes=True))

@pytest.fixture
def rerun_Gm(smallprojdir, run_small_T) -> RunDirs:
  """ Run the Gm mode of Zombi. """
  return RunDirs(run_small_T, G_runner(smallprojdir, mode='Gm', allgenomes=True))


@pytest.fixture(scope='session')
def run_Gu(projdir, run_T, run_RateCustomizer) -> RunDirs:
  """ Run the Gu mode of Zombi. Module-scoped to avoid rerunning. """
  assert run_RateCustomizer, 'There was a problem with run_RateCustomizer!'
  return RunDirs(run_T, G_runner(projdir, mode='Gu', allgenomes=True))

@pytest.fixture
def rerun_Gu(projdir, run_T, run_RateCustomizer) -> RunDirs:
  """ Run the Gu mode of Zombi. """
  assert run_RateCustomizer, 'There was a problem with run_RateCustomizer!'
  return RunDirs(run_T, G_runner(projdir, mode='Gu', allgenomes=True))

@pytest.fixture(scope='session')
def run_RateCustomizer(projdir):
  """ Test the RateCustomizer mode of Zombi. """
  result = subprocess.run(['zombiRateCustomizer', 'G', G_PARAMS_ALL, projdir])
  assert result.returncode == 0, (
    f'Error running RateCustomizer:\n{result.stderr}\n{result.stdout}')

  customrates = projdir / 'CustomRates'
  assert (customrates / TRANSFERRATES).exists()
  assert (customrates / EVENTRATES).exists()
  assert (customrates / EXTENSIONRATES).exists()
  return True
