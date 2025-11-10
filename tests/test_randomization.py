"""
Test that two runs of Zombi with the same seed, for each mode (T, G, S),
produce identical outputs.
"""
import filecmp
import pytest

from pathlib import Path

T_PARAMS = Path('tests/SpeciesTreeParametersSeeded.tsv')
G_PARAMS = Path('tests/GenomeParametersSeeded.tsv')
S_PARAMS = Path('tests/SequenceParametersSeeded.tsv')

NUMTHREADS = 28

@pytest.fixture(scope='session')
def basedir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project')


#@pytest.mark.dependency()
@pytest.fixture
def test_T(basedir, script_runner):
  """ Compare two runs of T mode using the same seed. """
  proj1 = basedir / 'proj1'
  script_runner.run(['zombi', 'T', T_PARAMS, proj1])

  proj2 = basedir / 'proj2'
  script_runner.run(['zombi', 'T', T_PARAMS, proj2])

  dcmp = filecmp.dircmp(proj1, proj2)
  assert identical_dirs(proj1 / 'T', proj2 / 'T'), ("T mode outputs differ: "
                                                    f"{dcmp.report_full_closure()}")

  return proj1 / 'T', proj2 / 'T'


#@pytest.mark.dependency(depends=['test_T'])
@pytest.fixture
def test_G(basedir, script_runner, test_T):
  """ Compare two runs of G mode using the same seed. """
  proj1 = basedir / 'proj1'
  assert (test_T[0]).exists(), 'Run test_T first!'
  script_runner.run(['zombi', 'G', G_PARAMS, proj1])

  proj2 = basedir / 'proj2'
  assert (test_T[1]).exists(), 'Run test_T first!'
  script_runner.run(['zombi', 'G', G_PARAMS, proj2])

  dcmp = filecmp.dircmp(proj1 / 'G', proj2 / 'G')
  assert identical_dirs(proj1 / 'G', proj2 / 'G'), ("G mode outputs differ: "
                                                    f"{dcmp.report_full_closure()}")

  return proj1 / 'G', proj2 / 'G'


@pytest.fixture
def test_Gf(basedir, script_runner, test_T):
  """ Compare two runs of G mode using the same seed. """
  proj1 = basedir / 'proj1'
  assert (test_T[0]).exists(), 'Run test_T first!'
  script_runner.run(['zombi', 'Gf', G_PARAMS, proj1])

  proj2 = basedir / 'proj2'
  assert (test_T[1]).exists(), 'Run test_T first!'
  script_runner.run(['zombi', 'Gf', G_PARAMS, proj2])

  dcmp = filecmp.dircmp(proj1 / 'G', proj2 / 'G')
  assert identical_dirs(proj1 / 'G', proj2 / 'G'), ("G mode outputs differ: "
                                                    f"{dcmp.report_full_closure()}")

  return proj1 / 'G', proj2 / 'G'


#@pytest.mark.dependency(depends=['test_T', 'test_G'])
def test_S(basedir, script_runner, test_G):
  """ Compare two runs of S mode using the same seed. """
  proj1 = basedir / 'proj1'
  assert (test_G[0]).exists(), 'Run test_G first!'
  script_runner.run(['zombi', 'S', f'-p {NUMTHREADS}', S_PARAMS, proj1])

  proj2 = basedir / 'proj2'
  assert (test_G[1]).exists(), 'Run test_G first!'
  script_runner.run(['zombi', 'S', f'-p {NUMTHREADS}', S_PARAMS, proj2])

  dcmp = filecmp.dircmp(proj1 / 'S', proj2 / 'S')
  assert identical_dirs(proj1 / 'S', proj2 / 'S'), ("S mode outputs differ: "
                                                    f"{dcmp.report_full_closure()}")


def test_Sf(basedir, script_runner, test_Gf):
  """ Compare two runs of S mode using the same seed. """
  proj1 = basedir / 'proj1'
  assert (test_Gf[0]).exists(), 'Run test_G first!'
  script_runner.run(['zombi', 'Sf', f'-p {NUMTHREADS}', S_PARAMS, proj1])

  proj2 = basedir / 'proj2'
  assert (test_Gf[1]).exists(), 'Run test_G first!'
  script_runner.run(['zombi', 'Sf', f'-p {NUMTHREADS}', S_PARAMS, proj2])

  dcmp = filecmp.dircmp(proj1 / 'S', proj2 / 'S')
  assert identical_dirs(proj1 / 'S', proj2 / 'S'), ("S mode outputs differ: "
                                                    f"{dcmp.report_full_closure()}")


#_______________________________________________________________________________
# Functions

def identical_dirs(dir1: Path, dir2: Path) -> bool:
  """
  Compare two directory trees content.
  Return False if they differ, True is they are the same.
  (https://stackoverflow.com/questions/4187564/recursively-compare-two-directories-to-ensure-they-have-the-same-files-and-subdi)
  """
  compare = filecmp.dircmp(dir1, dir2)
  if(compare.left_only or compare.right_only or compare.diff_files 
     or compare.funny_files):
    return False

  for subdir in compare.common_dirs:
    if not identical_dirs(dir1 / subdir, dir2 / subdir):
      return False

  return True
