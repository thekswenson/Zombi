"""
Unittests for testing the GenomeEvents and how they map coordinates.
"""

import os
import unittest
from zombi.Events import Inversion, LOSS, MapOriginError, MapPseudogeneError, ORIG, POS, TandemDup
from zombi.GenomeSimulator import GenomeSimulator
from zombi.Genomes import LEFT, RIGHT
import zombi.AuxiliarFunctions as af


GENOME_PARAMS = 'test/GenomeParametersDivisions.tsv'
TEST_FOLDER = 'test/TestDivisions/'
TEST_GENOME_30_6 = 'test/30_6.gff'  #30 bases, 5 * length-3 genomic/intergenomic pairs


class TestDivisions(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME_30_6):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER, 'T/Events.tsv')

    self.gss = GenomeSimulator(params, events_file, genome_file)
  
  def test_single_inversion(self):
    event1 = ("G", 1.3, "I", "n2", (2, 6, RIGHT))
    ## Events are a tuple where the elements are
    # 1. G or T (genome or tree level event)
    # 2. Time of the event
    # 3. Lineage undergoing the event
    # 4. Tuple with the details of the event
    
    self.gss.run_f_debug([event1])
    self.assertEqual(len(self.gss.all_genomes), 6) # We have six genomes


if __name__ == '__main__':
    unittest.main()
