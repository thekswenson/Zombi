"""
Unittests for testing the GenomeEvents and how they map coordinates.
"""

import os
import unittest
from zombi.Events import Inversion, TandemDup
from zombi.GenomeSimulator import GenomeSimulator
from zombi.Genomes import LEFT, RIGHT
import zombi.AuxiliarFunctions as af


GENOME_PARAMS = 'Parameters/GenomeParameters.tsv'
TEST_FOLDER = 'test/test_output/'
TEST_GENOME_30_10 = 'test/30_10.gff'  #30 bases, 3 * length-5 genomic/intergenomic pairs
TEST_GENOME_30_6 = 'test/30_6.gff'  #30 bases, 5 * length-3 genomic/intergenomic pairs
TEST_GENOME_18_6 = 'test/18_6.gff'  #18 bases, 3 * length-3 genomic/intergenomic pairs

class TestEvent(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME_30_6):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER, 'T/Events.tsv')

    self.gss = GenomeSimulator(params, events_file, genome_file)
    self.genome = self.gss.read_genome(genome_file, intergenic_sequences=True)

    self.gss.active_genomes.add(self.genome.species)
    self.gss.all_genomes["Root"] = self.genome

  def test_inversion_0(self):
    self.setUp(TEST_GENOME_30_10)

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
    self.setUp(TEST_GENOME_30_10)

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
    self.setUp(TEST_GENOME_30_10)

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

  def test_tandemdup_1(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_duplication_within_intergene(3, 9, RIGHT, lineage, 0.0)
    tdup: TandemDup = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 3,
                     'first intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[2].length, 1,
                     'second intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[4].length, 3,
                     'third intergene length mismatch after tandemdup')

    self.assertEqual(tdup.afterToBeforeS(3), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(8), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(9), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(10), 4,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(14), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(15), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(17), 11,
                     'intergene breakpoint mismap')

  def test_tandemdup_2RIGHT(self):
    self.setUp(TEST_GENOME_18_6)

    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_duplication_within_intergene(9, 2, RIGHT, lineage, 0.0)
    tdup: TandemDup = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 4,
                     'first intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[1].length, 3,
                     'second intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[3].length, 3,
                     'fourth intergene length mismatch after tandemdup')

    self.assertEqual(tdup.afterToBeforeS(2), 2,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(3), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(4), 11,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(5), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(6), 1,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(8), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(12), 7,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(14), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(16), 11,
                     'intergene breakpoint mismap')

    self.assertEqual(tdup.afterToBeforeT(0), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(5), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(6), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(7), 18,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(8), 1,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(12), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(23), 16,
                     'intergene breakpoint mismap')

  def test_tandemdup_2LEFT(self):
    self.setUp(TEST_GENOME_18_6)

    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_duplication_within_intergene(2, 9, LEFT, lineage, 0.0)
    tdup: TandemDup = ch.event_history[0]


    self.assertEqual(ch.intergenes[0].length, 4,
                     'first intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[1].length, 3,
                     'second intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[3].length, 3,
                     'fourth intergene length mismatch after tandemdup')

    self.assertEqual(tdup.afterToBeforeS(2), 2,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(3), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(4), 11,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(5), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(6), 1,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(8), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(12), 7,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(14), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(16), 11,
                     'intergene breakpoint mismap')

  def test_tandemdup_3(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_duplication_within_intergene(13, 8, RIGHT, lineage, 0.0)
    tdup: TandemDup = ch.event_history[0]


    self.assertEqual(ch.intergenes[2].length, 2,
                     'third intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[6].length, 3,
                     'seventh intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[7].length, 3,
                     'eigth intergene length mismatch after tandemdup')

    self.assertEqual(tdup.afterToBeforeS(8), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(9), 14,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(14), 19,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(15), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(22), 7,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(23), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(27), 12,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(30), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeS(34), 19,
                     'intergene breakpoint mismap')

    self.assertEqual(tdup.afterToBeforeT(0), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(7), 7,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(16), 23,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(23), 30,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(24), 1,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(38), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(42), 19,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(45), 22,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(48), 25,
                     'intergene breakpoint mismap')
    self.assertEqual(tdup.afterToBeforeT(53), 30,
                     'intergene breakpoint mismap')
 
if __name__ == '__main__':
    unittest.main()
