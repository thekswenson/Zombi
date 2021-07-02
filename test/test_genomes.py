"""
Unittests for testing the GenomeEvents and how they map coordinates.
"""

import os
import unittest
from zombi.GenomeSimulator import GenomeSimulator
from zombi.Genomes import LEFT, RIGHT
import zombi.AuxiliarFunctions as af


GENOME_PARAMS = 'Parameters/GenomeParameters.tsv'
TEST_FOLDER = 'test/test_output/'
TEST_GENOME_30_10 = 'test/100_10.gff'  #10 bases, 3 * length-5 genomic/intergenomic pairs
REPS = 1000

class TestGenomes(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME_30_10):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER, 'T/Events.tsv')

    self.gss = GenomeSimulator(params, events_file, genome_file)
    self.genome = self.gss.read_genome(genome_file, intergenic_sequences=True)

    self.gss.active_genomes.add(self.genome.species)
    self.gss.all_genomes["Root"] = self.genome

  def test_coordinate_selection_1(self):
    ch = self.genome.chromosomes[0]

    chosen = []
    for _ in range(REPS):
        # The specific coordinate ranges for the intergenes are the following:
        # 0-5, 6-11, 12-17, 18-23, 24-29, 30-35, 36-41, 42-47, 48-53, 54-59
      c = ch.select_random_intergenic_coordinate_excluding(16, 38, RIGHT)
      self.assertTrue(0 <= c <= 11 or 42 <= c <= 59, f'bad coordinate: {c}')
      chosen.append(c)

    self.assertTrue(set(chosen) == set(range(12)) | set(range(42, 60)),
                    f'low probability event occured (coupon collectors problem)')

    

if __name__ == '__main__':
    unittest.main()
