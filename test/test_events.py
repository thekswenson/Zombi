"""
Unittests for testing the GenomeEvents and how they map coordinates.
"""

import os
import unittest
from zombi.Events import Inversion, Loss, Origination, TandemDup, Transposition
from zombi.Events import MapOriginError, MapPseudogeneError, Transfer
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

  # INVERSIONS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def test_inversion_0(self):
    self.setUp(TEST_GENOME_30_10)

    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_inversion_intergenic(ch, 0, 6, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 0,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[0].length, inversion.afterL.specificLen()-1,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[1].length, 10,
                     'second intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[1].length, inversion.afterR.specificLen()-1,
                     'second intergene length mismatch after inversion')

      #Specific coordinate mapping:
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

      #Total coordinate mapping:
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
    self.gss.make_inversion_intergenic(ch, 3, 10, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 7,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[1].length, 3,
                     'second intergene length mismatch after inversion')

      #Specific coordinate mapping:
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

      #Total coordinate mapping:
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
    self.gss.make_inversion_intergenic(ch, 3, 10, LEFT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 3,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[1].length, 7,
                     'second intergene length mismatch after inversion')

      #Specific coordinate mapping:
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

      #Total coordinate mapping:
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
    self.gss.make_inversion_intergenic(ch, 9, 5, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[1].length, 4,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[2].length, 2,
                     'second intergene length mismatch after inversion')

      #Specific coordinate mapping:
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

      #Total coordinate mapping:
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
    self.gss.make_inversion_intergenic(ch, 13, 8, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[2].length, 5,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[3].length, 1,
                     'second intergene length mismatch after inversion')

      #Specific coordinate mapping:
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

      #Total coordinate mapping:
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

  def test_inversion_Adri_19_10_21(self):
    """
    This test wraps.
    """
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_inversion_intergenic(ch, 16, 10, RIGHT, lineage, 0.0)
    inversion: Inversion = ch.event_history[0]

    self.assertEqual(ch.intergenes[2].length, 4,
                     'first intergene length mismatch after inversion')
    self.assertEqual(ch.intergenes[4].length, 2,
                     'second intergene length mismatch after inversion')

      #Specific coordinate mapping:
    self.assertEqual(inversion.afterToBeforeS(0), 7,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(10), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(14), 13,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(19), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(18), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeS(17), 16,
                     'intergene breakpoint mismap')

      #Total coordinate mapping:
    self.assertEqual(inversion.afterToBeforeT(3), 12,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(17), 28,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(23), 22,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(30), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(29), 16,
                     'intergene breakpoint mismap')
    self.assertEqual(inversion.afterToBeforeT(28), 27,
                     'intergene breakpoint mismap')

  # TANDEM DUPLICATIONS - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def test_tandemdup_1(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_duplication_within_intergene(ch, 3, 9, RIGHT, lineage, 0.0)
    tdup: TandemDup = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 3,
                     'first intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[2].length, 1,
                     'second intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[4].length, 3,
                     'third intergene length mismatch after tandemdup')

      #Specific coordinate mapping:
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
    self.gss.make_duplication_within_intergene(ch, 9, 2, RIGHT, lineage, 0.0)
    tdup: TandemDup = ch.event_history[0]

    self.assertEqual(ch.intergenes[0].length, 4,
                     'first intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[1].length, 3,
                     'second intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[3].length, 3,
                     'fourth intergene length mismatch after tandemdup')

      #Specific coordinate mapping:
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

      #Total coordinate mapping:
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
    self.gss.make_duplication_within_intergene(ch, 2, 9, LEFT, lineage, 0.0)
    tdup: TandemDup = ch.event_history[0]


    self.assertEqual(ch.intergenes[0].length, 4,
                     'first intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[1].length, 3,
                     'second intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[3].length, 3,
                     'fourth intergene length mismatch after tandemdup')

      #Specific coordinate mapping:
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
    self.gss.make_duplication_within_intergene(ch, 13, 8, RIGHT, lineage, 0.0)
    tdup: TandemDup = ch.event_history[0]


    self.assertEqual(ch.intergenes[2].length, 2,
                     'third intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[6].length, 3,
                     'seventh intergene length mismatch after tandemdup')
    self.assertEqual(ch.intergenes[7].length, 3,
                     'eigth intergene length mismatch after tandemdup')

      #Specific coordinate mapping:
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

      #Total coordinate mapping:
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

  # LOSSES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def test_loss_1(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_loss_intergenic(ch, 3, 10, RIGHT, lineage, 0.0)
    loss: Loss = ch.event_history[0]


    self.assertEqual(ch.intergenes[0].length, 4,
                     'first intergene length mismatch after loss')

      #Specific coordinate mapping:
    self.assertEqual(loss.afterToBeforeS(3), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(4), 11,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(12), 19,
                     'intergene breakpoint mismap')

      #Total coordinate mapping:
    self.assertEqual(loss.afterToBeforeT(6), 6,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(7), 18,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(15), 26,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(19), 30,
                     'intergene breakpoint mismap')

  def test_loss_1P(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_loss_intergenic(ch, 3, 10, RIGHT, lineage, 0.0, True)
    loss: Loss = ch.event_history[0]


    self.assertEqual(ch.intergenes[0].length, 15,
                     'first intergene length mismatch after psuedogenization')

      #Specific coordinate mapping:
    self.assertEqual(loss.afterToBeforeS(3), 3,
                     'intergene breakpoint mismap')
    with self.assertRaises(MapPseudogeneError):
      loss.afterToBeforeS(4)
    self.assertEqual(loss.afterToBeforeS(6), 4,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(9), 7,
                     'intergene breakpoint mismap')
    with self.assertRaises(MapPseudogeneError):
      loss.afterToBeforeS(10)
    with self.assertRaises(MapPseudogeneError):
      loss.afterToBeforeS(11)
    self.assertEqual(loss.afterToBeforeS(12), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(13), 9,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(14), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(23), 19,
                     'intergene breakpoint mismap')

      #Total coordinate mapping:
    self.assertEqual(loss.afterToBeforeT(6), 6,
                     'intergene breakpoint mismap')
    with self.assertRaises(NotImplementedError):
      loss.afterToBeforeT(7)
    with self.assertRaises(NotImplementedError):
      loss.afterToBeforeT(16)
    self.assertEqual(loss.afterToBeforeT(17), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(17), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(30), 30,
                     'intergene breakpoint mismap')

  def test_loss_2(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_loss_intergenic(ch, 14, 6, RIGHT, lineage, 0.0)
    loss: Loss = ch.event_history[0]


    self.assertEqual(ch.intergenes[1].length, 3,
                     'second intergene length mismatch after loss')

    self.assertEqual(loss.afterToBeforeS(0), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(6), 14,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(7), 7,
                     'intergene breakpoint mismap')

    self.assertEqual(loss.afterToBeforeT(0), 12,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(3), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(11), 23,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(12), 12,
                     'intergene breakpoint mismap')

  def test_loss_2P(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_loss_intergenic(ch, 14, 6, RIGHT, lineage, 0.0, True)
    loss: Loss = ch.event_history[0]


    self.assertEqual(ch.intergenes[1].length, 21,
                     'second intergene length mismatch after pseudogenization')

      #Specific coordinate mapping:
    self.assertEqual(loss.afterToBeforeS(0), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(5), 13,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(6), 14,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(7), 15,
                     'intergene breakpoint mismap')
    with self.assertRaises(MapPseudogeneError):
      loss.afterToBeforeS(8)
    with self.assertRaises(MapPseudogeneError):
      loss.afterToBeforeS(9)
    self.assertEqual(loss.afterToBeforeS(10), 16,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(13), 19,
                     'intergene breakpoint mismap')
    with self.assertRaises(MapPseudogeneError):
      loss.afterToBeforeS(14)
    self.assertEqual(loss.afterToBeforeS(16), 0,
                     'intergene breakpoint mismap')
    with self.assertRaises(MapPseudogeneError):
      loss.afterToBeforeS(21)
    self.assertEqual(loss.afterToBeforeS(22), 4,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeS(25), 7,
                     'intergene breakpoint mismap')

      #Total coordinate mapping:
    self.assertEqual(loss.afterToBeforeT(0), 12,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(11), 23,
                     'intergene breakpoint mismap')
    with self.assertRaises(NotImplementedError):
      loss.afterToBeforeT(12)
    with self.assertRaises(NotImplementedError):
      loss.afterToBeforeT(28)
    self.assertEqual(loss.afterToBeforeT(29), 11,
                     'intergene breakpoint mismap')
    self.assertEqual(loss.afterToBeforeT(30), 12,
                     'intergene breakpoint mismap')

  # TRANSPOSITIONS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def test_transposition1(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_transposition_intergenic(ch, 2, 9, RIGHT, 17, lineage, 0.0)
    trans: Transposition = ch.event_history[0]


    self.assertEqual(ch.intergenes[0].length, 4,
                     'first intergene length mismatch after loss')
    self.assertEqual(ch.intergenes[2].length, 2,
                     'third intergene length mismatch after loss')
    self.assertEqual(ch.intergenes[4].length, 3,
                     'fifth intergene length mismatch after loss')

      #Specific coordinate mapping:
    self.assertEqual(trans.afterToBeforeS(2), 2,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(3), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(10), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(11), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(14), 6,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(16), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(17), 17,
                     'intergene breakpoint mismap')

      #Total coordinate mapping:
    self.assertEqual(trans.afterToBeforeT(5), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(6), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(17), 28,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(18), 6,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(27), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(28), 28,
                     'intergene breakpoint mismap')

  def test_transposition2(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_transposition_intergenic(ch, 6, 11, RIGHT, 0, lineage, 0.0)
    trans: Transposition = ch.event_history[0]


    self.assertEqual(ch.intergenes[0].length, 1,
                     'first intergene length mismatch after loss')
    self.assertEqual(ch.intergenes[1].length, 6,
                     'second intergene length mismatch after loss')
    self.assertEqual(ch.intergenes[2].length, 2,
                     'third intergene length mismatch after loss')

      #Specific coordinate mapping:
    self.assertEqual(trans.afterToBeforeS(0), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(1), 7,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(4), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(5), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(10), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(11), 11,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(14), 14,
                     'intergene breakpoint mismap')

      #Total coordinate mapping:
    self.assertEqual(trans.afterToBeforeT(3), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(4), 12,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(9), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(10), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(15), 8,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(17), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(18), 18,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(25), 25,
                     'intergene breakpoint mismap')

  def test_transposition3(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    self.gss.make_transposition_intergenic(ch, 14, 1, RIGHT, 10, lineage, 0.0)
    trans: Origination = ch.event_history[0]


    self.assertEqual(ch.intergenes[0].length, 3,
                     'first intergene length mismatch after loss')
    self.assertEqual(ch.intergenes[1].length, 3,
                     'second intergene length mismatch after loss')
    self.assertEqual(ch.intergenes[3].length, 2,
                     'fourth intergene length mismatch after loss')
    self.assertEqual(ch.intergenes[4].length, 4,
                     'fifth intergene length mismatch after loss')

      #Specific coordinate mapping:
    self.assertEqual(trans.afterToBeforeS(0), 4,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(6), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(7), 15,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(11), 19,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(12), 0,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(13), 10,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(17), 14,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS(18), 2,
                     'intergene breakpoint mismap')

      #Total coordinate mapping:
    self.assertEqual(trans.afterToBeforeT(0), 6,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(11), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(12), 24,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(18), 30,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(19), 1,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(21), 3,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(22), 17,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(28), 23,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(29), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT(30), 6,
                     'intergene breakpoint mismap')

  # ORIGINATIONS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def test_origination1(self):
    ch = self.genome.chromosomes[0]
    lineage = self.genome.species
    gene = self.gss.make_origination_intergenic(ch, 2, lineage, 0.0)
    orig: Origination = ch.event_history[0]


    self.assertEqual(ch.intergenes[0].length, 2,
                     'first intergene length mismatch after loss')
    self.assertEqual(ch.intergenes[1].length, 1,
                     'second intergene length mismatch after loss')

      #Specific coordinate mapping:
    self.assertEqual(orig.afterToBeforeS(2), 2,
                     'intergene breakpoint mismap')
    self.assertEqual(orig.afterToBeforeS(3), 2,
                     'intergene breakpoint mismap')
    self.assertEqual(orig.afterToBeforeS(4), 3,
                     'intergene breakpoint mismap')

      #Total coordinate mapping:
    self.assertEqual(orig.afterToBeforeT(5), 5,
                     'intergene breakpoint mismap')
    with self.assertRaises(MapOriginError):
      orig.afterToBeforeT(6)
    with self.assertRaises(MapOriginError):
      orig.afterToBeforeT(4 + gene.length)
    self.assertEqual(orig.afterToBeforeT(5 + gene.length), 5,
                     'intergene breakpoint mismap')
    self.assertEqual(orig.afterToBeforeT(6 + gene.length), 6,
                     'intergene breakpoint mismap')

  # TRANSFERS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def transferSetUp(self):
    donorlineage = 'donorlineage'
    receptorlineage = 'receptorlineage'
    g1, g2 = self.gss.make_speciation(self.genome.species, donorlineage,
                                      receptorlineage, 1.0)
    donor = g1.chromosomes[0]
    receptor = g2.chromosomes[0]

    donor.obtain_locations()
    receptor.obtain_locations()

    return donor, donorlineage, receptor, receptorlineage

  def test_transfer1(self):
    self.setUp(TEST_GENOME_18_6)
    donor, donorlineage, receptor, receptorlineage = self.transferSetUp()

    self.gss.make_transfer_intergenic(donor, 2, 10, RIGHT, donorlineage,
                                      receptor, 1, receptorlineage, 1.5)
    trans: Transfer = receptor.event_history[0]

    self.assertEqual(trans.afterToBeforeS_lineage(1), (receptorlineage, 1),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(2), (donorlineage, 3),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(8), (donorlineage, 9),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(9), (receptorlineage, 1),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(10), (receptorlineage, 2),
                     'intergene breakpoint mismap')

    self.assertEqual(trans.afterToBeforeT_lineage(4), (receptorlineage, 4),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(5), (donorlineage, 6),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(12), (donorlineage, 13),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(15), (donorlineage, 16),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(16), (receptorlineage, 4),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(17), (receptorlineage, 5),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(26), (receptorlineage, 14),
                     'intergene breakpoint mismap')

  def test_transfer2(self):
    self.setUp(TEST_GENOME_18_6)
    donor, donorlineage, receptor, receptorlineage = self.transferSetUp()

    self.gss.make_transfer_intergenic(donor, 3, 10, LEFT, donorlineage,
                                      receptor, 7, receptorlineage, 1.5)
    trans: Transfer = receptor.event_history[0]

    self.assertEqual(trans.afterToBeforeS_lineage(7), (receptorlineage, 7),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(8), (donorlineage, 11),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(9), (donorlineage, 0),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(11), (donorlineage, 2),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(12), (receptorlineage, 7),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeS_lineage(15), (receptorlineage, 10),
                     'intergene breakpoint mismap')

    self.assertEqual(trans.afterToBeforeT_lineage(12), (receptorlineage, 12),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(13), (donorlineage, 18),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(14), (donorlineage, 1),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(18), (donorlineage, 5),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(19), (receptorlineage, 12),
                     'intergene breakpoint mismap')
    self.assertEqual(trans.afterToBeforeT_lineage(22), (receptorlineage, 15),
                     'intergene breakpoint mismap')

if __name__ == '__main__':
    unittest.main()
