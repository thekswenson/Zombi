"""
Unittests for testing the GenomeEvents and how they map coordinates.
"""

import os
import unittest
from zombi.Events import Inversion
from zombi.GenomeSimulator import GenomeSimulator
from zombi.Genomes import LEFT, RIGHT
import zombi.AuxiliarFunctions as af


GENOME_PARAMS = 'Parameters/GenomeParameters.tsv'
TEST_FOLDER = 'test/small_output/'
TEST_GENOME_30_10 = 'test/30_10.gff'  #30 bases, 3 * length-5 genomic/intergenomic pairs
TEST_GENOME_30_6 = 'test/30_6.gff'  #30 bases, 5 * length-3 genomic/intergenomic pairs

class TestEvent(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME_30_10):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER, 'T/Events.tsv')

    self.gss = GenomeSimulator(params, events_file, genome_file)
    self.genome = self.gss.read_genome(genome_file, intergenic_sequences=True)

    self.gss.active_genomes.add(self.genome.species)
    self.gss.all_genomes["Root"] = self.genome

  def test_inversion_0(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_inversion_intergenic(0, 6, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 0,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[0].length, inversion.afterL.specificLen()-1,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[1].length, 10,
                     'second intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[1].length, inversion.afterR.specificLen()-1,
                     'second intergene length mismatch after inversion')

    self.assertEqual(inversion.afterToBeforeS(0), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(1), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(4), 2,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(6), 6,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(10), 10,
                     'intergene breakpoint mismap')

    self.assertEqual(inversion.afterToBeforeT(0), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(5), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(6), 14,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(10), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(14), 6,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(15), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(20), 20,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(30), 30,
                     'intergene breakpoint mismap')

  def test_inversion_1RIGHT(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_inversion_intergenic(3, 10, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 7,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[1].length, 3,
                     'second intergene length mismatch after inversion')

    self.assertEqual(inversion.afterToBeforeS(0), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(3), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(6), 7,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(8), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(10), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(13), 13,
                     'intergene breakpoint mismap')

    self.assertEqual(inversion.afterToBeforeT(8), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(9), 18,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(15), 12,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(18), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(19), 19,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(22), 22,
                     'intergene breakpoint mismap')

  def test_inversion_1LEFT(self):
    """
    This test wraps.
    """
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_inversion_intergenic(3, 10, LEFT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 3,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[1].length, 7,
                     'second intergene length mismatch after inversion')

    self.assertEqual(inversion.afterToBeforeS(0), 11,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(1), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(6), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(8), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(9), 2,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(11), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(12), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(17), 12,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(14), 15,
                     'intergene breakpoint mismap')

    self.assertEqual(inversion.afterToBeforeT(0), 25,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(5), 20,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(6), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(9), 11,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(17), 19,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(18), 7,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(24), 1,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(25), 0,
                     'intergene breakpoint mismap')

  def test_inversion_2(self):
    """
    This test wraps.
    """
    self.setUp(TEST_GENOME_30_6)

    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_inversion_intergenic(9, 5, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[1].length, 4,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[2].length, 2,
                     'second intergene length mismatch after inversion')

    self.assertEqual(inversion.afterToBeforeS(0), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(5), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(6), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(10), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(11), 4,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(12), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(15), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(16), 19,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(19), 16,
                     'intergene breakpoint mismap')

    self.assertEqual(inversion.afterToBeforeT(0), 27,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(5), 22,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(11), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(14), 13,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(17), 16,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(18), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(25), 2,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(29), 28,
                     'intergene breakpoint mismap')

  def test_inversion_3(self):
    """
    This test wraps.
    """
    self.setUp(TEST_GENOME_30_6)

    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_inversion_intergenic(13, 8, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[2].length, 5,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[3].length, 1,
                     'second intergene length mismatch after inversion')

    self.assertEqual(inversion.afterToBeforeS(0), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(4), 19,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(7), 16,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(9), 14,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(10), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(15), 13,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(16), 7,
                     'intergene breakpoint mismap')

    self.assertEqual(inversion.afterToBeforeT(0), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(6), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(8), 1,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(9), 30,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(10), 29,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(16), 23,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(17), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(21), 19,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(24), 22,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(25), 14,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(30), 9,
                     'intergene breakpoint mismap')

if __name__ == '__main__':
    unittest.main()
