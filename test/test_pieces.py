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
TEST_GENOME = 'test/30_6.gff'


class TestDivisions(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER1, 'T/Events.tsv')

    self.gss = GenomeSimulator(params, events_file, genome_file)
  
  def test_divisions(self):

      event1 = ("G", 1.22, "I", "n1", (16, 10, RIGHT))
      event2 = ("G", 1.233, "I", "n1", (3, 14, RIGHT))

      #self.gss.run_f()
      #self.gss.run_f_debug([event1,])
      self.gss.run_f_debug([event1, event2])  
      self.gss.obtain_divisions()       
      self.gss.obtain_events_for_divisions()

      #for ch in self.gss.all_genomes_second["n1"]:
      #    ch.print_pieces()      

      
      
      print("Should look like:")
      for ch in self.gss.all_genomes["n1"]:
          for gene, intergene in zip(ch.genes, ch.intergenes):
              print("Gene", gene.total_flanking, gene.gene_family)
              print("Intergene", intergene.total_flanking)

      #print("**")

      
      

if __name__ == '__main__':
    unittest.main()
