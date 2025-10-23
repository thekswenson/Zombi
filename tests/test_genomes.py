"""
Unittests for testing the GenomeEvents and how they map coordinates.
"""
import shutil
import pytest
import unittest
from pathlib import Path

from zombi.Filenames import TREEEVENTS
from zombi.GenomeSimulator import GenomeSimulator
from zombi.Genomes import RIGHT
import zombi.AuxiliarFunctions as af


GENOME_PARAMS = Path('Parameters/GenomeParameters.tsv')
TEST_DIVISIONS = Path('tests/TestDivisions1/')
TEST_GENOME_30_10 = Path('tests/100_10.gff')  #10 bases, 3 * length-5 genomic/intergenomic pairs
REPS = 1000

@pytest.fixture(autouse=True)
def _inject_tmp_path_factory(request, tmp_path_factory):
    if hasattr(request, "cls") and request.cls is not None:
        request.cls.tmp_path_factory = tmp_path_factory


@pytest.mark.usefixtures("tmp_path_factory")
class TestGenomes(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME_30_10):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    test_folder = self.tmp_path_factory.mktemp('session')   #type: ignore

    shutil.copytree(TEST_DIVISIONS / 'T', test_folder / 'T')
    events_file = test_folder / 'T' / TREEEVENTS

    self.gss = GenomeSimulator(params, events_file, genome_file)
    self.genome = self.gss.read_genome(genome_file, intergenic_sequences=True)

    self.gss.active_genomes.add(self.genome.species)
    self.gss.node_genomes["Root"] = self.genome

  def test_coordinate_selection_1(self):
    ch = self.genome.chromosomes[0]

    chosen = []
    for _ in range(REPS):
        # The specific coordinate ranges for the intergenes are the following:
        # 0-5, 6-11, 12-17, 18-23, 24-29, 30-35, 36-41, 42-47, 48-53, 54-59
      c = ch.select_random_intergenic_coordinate_excluding(16, 38, RIGHT)
      assert c is not None
      self.assertTrue(0 <= c <= 11 or 42 <= c <= 59, f'bad coordinate: {c}')
      chosen.append(c)

    self.assertTrue(set(chosen) == set(range(12)) | set(range(42, 60)),
                    f'low probability event occured (coupon collectors problem)')

    

if __name__ == '__main__':
    unittest.main()
