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

      event0 = ("G", 1.21, "I", "n1", (1, 9, RIGHT))
      event1 = ("G", 1.22, "I", "n1", (15, 6, RIGHT))
      event2 = ("G", 1.233, "I", "n1", (11, 1, RIGHT))
      event3 = ("G", 1.234, "I", "n1", (8, 17, RIGHT))

      #
      #self.gss.run_f()
      ##self.gss.run_f_debug([event1,event2])
      #self.gss.run_f_debug([event0, event1, event2, event3])  
      #self.gss.obtain_divisions()       
      #self.gss.obtain_events_for_divisions()

      #for ch in self.gss.all_genomes_second["n1"]:
      #    ch.print_pieces()      

      
      
      #print("Should look like:")
      #for ch in self.gss.all_genomes["n1"]:
      #    for gene, intergene in zip(ch.genes, ch.intergenes):
      #        print("Gene", gene.total_flanking, gene.gene_family)
      #        print("Intergene", intergene.total_flanking)

      #print("**")

  def test_divisions2(self):

      event0 = ("G", 0.028, "I", "Root", (5, 9, RIGHT))
      event1 = ("G", 0.21, "I", "Root", (19, 1, RIGHT))
      event2 = ("G", 0.37, "I", "Root", (8, 12, RIGHT))
      event3 = ("G", 0.40, "I", "Root", (10, 6, RIGHT))
      event4 = ("G", 0.44, "I", "Root", (19, 4, RIGHT))
      event5 = ("G", 0.48, "I", "Root", (5, 17, RIGHT))
      event6 = ("G", 0.53, "I", "Root", (18, 14, RIGHT))
      event7 = ("G", 0.644, "I", "Root", (11, 2, RIGHT))
      event8 = ("G", 0.74, "I", "Root", (14, 2, RIGHT))
      event9 = ("G", 0.748, "I", "Root", (12, 15, RIGHT))
      event10 = ("G", 0.919, "I", "Root", (2, 3, RIGHT))
      event11 = ("G", 1.27, "I", "n2", (1, 14, RIGHT))
      event12 = ("G", 1.42, "I", "n2", (13, 11, RIGHT))
      event13 = ("G", 1.48, "I", "n2", (6, 19, RIGHT))
      event14 = ("G", 1.55, "I", "n2", (0, 18, RIGHT))
      event15 = ("G", 1.47, "I", "n1", (13, 18, RIGHT))

      
      self.gss.run_f_debug([event1,event2])
      #self.gss.run_f_debug([event1,  event4, event9])  
      self.gss.obtain_divisions()       
      self.gss.obtain_events_for_divisions()

      #for ch in self.gss.all_genomes_second["n1"]:
      #    ch.print_pieces()      
   
  def test_transpositions1(self):
     
      event1 = ("G", 0.028, "P", "Root", (5, 9, RIGHT))
      #self.gss.run_f() 
      #self.gss.obtain_divisions()       
      #self.gss.obtain_events_for_divisions()

  def test_genetrees(self):
     
      event1 = ("G", 0.028, "P", "Root", (5, 9, RIGHT))
      
      self.gss.run_f() 
      self.gss.obtain_divisions()       
      self.gss.obtain_events_for_divisions()
      
      
      division_trees_folder = os.path.join('test/TestDivisions1/', "Division_trees")
      
      self.gss.write_division_trees('test/TestDivisions1/')
      self.gss.write_division_coordinates('test/TestDivisions1/')
      print(self.gss.all_division_families["1"].events)

if __name__ == '__main__':
    unittest.main()
