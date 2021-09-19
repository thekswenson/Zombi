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
TEST_FOLDER2 = 'test/TestDivisions2/'
TEST_GENOME_30_6 = 'test/30_6.gff'  #30 bases, 5 * length-3 genomic/intergenomic pairs


class TestDivisions1(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME_30_6):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER1, 'T/Events.tsv')

    self.gss = GenomeSimulator(params, events_file, genome_file)
  
  def test_divisions(self):
      self.gss.run_f_debug([])
      self.gss.obtain_divisions() # This obtain the divisions at the root
      self.assertEqual(self.gss.natural_cuts,
                       [(0,3),(4,7),(8,11),(12,15),(16,19)])

  def test_single_inversion_RIGHT(self):
    event1 = ("G", 1.3, "I", "n2", (2, 6, RIGHT))
    ## Events are a tuple where the elements are
    # 1. G or T (genome or tree level event)
    # 2. Time of the event
    # 3. Lineage undergoing the event
    # 4. Tuple with the details of the event
    
    self.gss.run_f_debug([event1])
    self.gss.obtain_divisions() # This obtain the divisions at the root
    self.assertEqual(self.gss.initial_divisions,
                     [(0,2),(2,3),(4,6),(6,7),(8,11),(12,15),(16,19)])
    
  def test_single_inversion_LEFT(self):
    event1 = ("G", 1.3, "I", "n2", (6, 2, LEFT))
    
    self.gss.run_f_debug([event1])
    self.gss.obtain_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,2),(2,3),(4,6),(6,7),(8,11),(12,15),(16,19)])

  def test_single_inversion_LEFT_wrapping(self): 
    event1 = ("G", 1.3, "I", "n2", (2, 6, LEFT))
    
    self.gss.run_f_debug([event1])
    self.gss.obtain_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,2),(2,3),(4,6),(6,7),(8,11),(12,15),(16,19)])

  def test_2inversions_1branch(self):
    event1 = ("G", 1.3, "I", "n2", (2, 6, LEFT))
    event2 = ("G", 1.31, "I", "n2", (4, 10, RIGHT))
    
    self.gss.run_f_debug([event1, event2])
    self.gss.obtain_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,2),(2,3),(4,5),(5,6),(6,7),(8,11),(12,15),(16,17),
                      (17,19)])
    
  def test_3inversions_1branch(self):
    event1 = ("G", 1.3, "I", "n2", (2, 6, LEFT))
    event2 = ("G", 1.31, "I", "n2", (4, 10, RIGHT))
    event3 = ("G", 1.32, "I", "n2", (8, 17, RIGHT))
    
    self.gss.run_f_debug([event1, event2, event3])
    self.gss.obtain_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,1),(1,2),(2,3),(4,5),(5,6),(6,7),(8,10),(10,11),
                      (12,15),(16,17),(17,19)])
  

class TestDivisions2(unittest.TestCase): # In a slightly more compex tree

  def setUp(self, genome_file=TEST_GENOME_30_6):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER2, 'T/Events.tsv')
    self.gss = GenomeSimulator(params, events_file, genome_file)

  def test_inversions_multiple_branches(self):

    event1 = ("G", 0.624, "I", "Root", (19, 1, RIGHT))
    event2 = ("G", 1.073, "I", "n1", (15, 5, LEFT))

    self.gss.run_f_debug([event1, event2])
    self.gss.obtain_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,1),(1,3),(4,6),(6,7),(8,11),(12,15),(16,19)])
  
  def test_events_for_divisions(self):

    print("***")

    #events_file = "/Users/aadavin/Desktop/Krister/Zombi/test/Events.tsv"
    
    #events = self.gss.read_genome_events_file(events_file)

    event1 = ("G", 0.024, "D", "Root", (2,12, RIGHT))
    event2 = ("G", 0.025, "I", "Root", (2, 6, RIGHT))
    event3 = ("G", 0.026, "D", "Root", (3, 7, RIGHT))
    event4 = ("G", 0.027, "I", "Root", (2, 11, RIGHT))
    event5 = ("G", 0.028, "D", "Root", (9, 15, RIGHT))
    event6 = ("G", 0.029, "I", "Root", (8, 15, RIGHT))

  
    #self.gss.run_f_debug([event1, event2, event3, event4, event5, event6])  

    #self.gss.run_f_debug([event1]) 
    # 
    #  

    self.gss.run_f()
    self.gss.obtain_divisions() 
    self.gss.obtain_events_for_divisions() 
  
    for ch in self.gss.all_genomes_second["n1"]:
         for intergene in ch.iter_intergenes():
             if len(intergene) == 0:
               print(intergene)
             for division in intergene:
                 print(intergene, len(intergene), division, len(division))

    #self.assertEqual()

if __name__ == '__main__':
    unittest.main()
