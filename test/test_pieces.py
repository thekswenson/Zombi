"""
Unittests for testing the GenomeEvents and how they map coordinates.
"""

import os
import unittest
from zombi.GenomeSimulator import GenomeSimulator
from zombi.Genomes import LEFT, RIGHT
import zombi.AuxiliarFunctions as af


GENOME_PARAMS = 'test/GenomeParametersDivisions.tsv'
TEST_FOLDER1 = 'test/TestDivisions1/'
TEST_GENOME_18_6 = 'test/18_6.gff'  #30 bases, 5 * length-3 genomic/intergenomic pairs


class TestDivisions(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME_18_6):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER1, 'T/Events.tsv')

    self.gss = GenomeSimulator(params, events_file, genome_file)
  
  def test_divisions(self):

      #event1 = ("G", 1.22, "I", "n1", (0, 11, RIGHT))
      #event1 = ("G", 1.23, "I", "n1", (11, 1, RIGHT))
      
      #event3 = ("G", 1.25, "I", "n1", (9, 4, RIGHT))
      #event4 = ("G", 1.26, "I", "n1", (11, 2, RIGHT))
      
      self.gss.run_f()

      #self.gss.run_f_debug([event1])
      self.gss.obtain_divisions() 
      self.gss.obtain_events_for_divisions()

      for ch in self.gss.all_genomes["n1"]:
          for gene, intergene in zip(ch.genes, ch.intergenes):
              print("Gene", gene.total_flanking, gene.gene_family)
              print("Intergene", intergene.total_flanking)

      
      print("***")
      for ch in self.gss.all_genomes_second["n1"]:
          ch.print_pieces()
      
      print("**")

      
      

if __name__ == '__main__':
    unittest.main()
