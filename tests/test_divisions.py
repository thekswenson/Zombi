"""
Pytests for testing the GenomeEvents and how they map coordinates.
"""
import unittest
from pathlib import Path

from zombi.Events import INV, TDUP
from zombi.GenomeSimulator import GenomeSimulator
from zombi.Genomes import T_DIR
import zombi.AuxiliarFunctions as af


GENOME_PARAMS = Path('tests/GenomeParametersDivisions.tsv')
TEST_GENOME_30_6 = Path('tests/30_6.gff')  #30 bases, 5 * length-3 genomic/intergenomic pairs
TEST_FOLDER1 = Path('tests/TestDivisions1')
TEST_FOLDER2 = Path('tests/TestDivisions2')


class TestDivisions1(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME_30_6):
    params = af.prepare_genome_parameters(GENOME_PARAMS)
    #events_file = self.projbase / TEST_FOLDER1 / 'T/Events.tsv'
    events_file = TEST_FOLDER1 / 'T/Events.tsv'

    self.gss = GenomeSimulator(params, events_file, genome_file)
  
  def test_divisions(self):
      self.gss.run_f_debug([])
      self.gss.init_divisions() # This obtain the divisions at the root
      self.assertEqual(self.gss.natural_cuts,
                       [(0,3),(4,7),(8,11),(12,15),(16,19)])

  def test_single_inversion_RIGHT(self):
    event1 = ("G", 1.3, INV, "n2", (2, 6, T_DIR.RIGHT))
    ## Events are a tuple where the elements are
    # 1. G or T (genome or tree level event)
    # 2. Time of the event
    # 3. Lineage undergoing the event
    # 4. Tuple with the details of the event
    
    self.gss.run_f_debug([event1])
    self.gss.init_divisions() # This obtain the divisions at the root
    self.assertEqual(self.gss.initial_divisions,
                     [(0,2),(2,3),(4,6),(6,7),(8,11),(12,15),(16,19)])
    
  def test_single_inversion_LEFT(self):
    event1 = ("G", 1.3, INV, "n2", (6, 2, T_DIR.LEFT))

    self.gss.run_f_debug([event1])
    self.gss.init_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,2),(2,3),(4,6),(6,7),(8,11),(12,15),(16,19)])

  def test_single_inversion_LEFT_wrapping(self): 
    event1 = ("G", 1.3, INV, "n2", (2, 6, T_DIR.LEFT))

    self.gss.run_f_debug([event1])
    self.gss.init_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,2),(2,3),(4,6),(6,7),(8,11),(12,15),(16,19)])

  def test_2inversions_1branch(self):
    event1 = ("G", 1.3, INV, "n2", (2, 6, T_DIR.LEFT))
    event2 = ("G", 1.31, INV, "n2", (4, 10, T_DIR.RIGHT))

    self.gss.run_f_debug([event1, event2])
    self.gss.init_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,2),(2,3),(4,5),(5,6),(6,7),(8,11),(12,15),(16,17),
                      (17,19)])
    
  def test_3inversions_1branch(self):
    event1 = ("G", 1.3, INV, "n2", (2, 6, T_DIR.LEFT))
    event2 = ("G", 1.31, INV, "n2", (4, 10, T_DIR.RIGHT))
    event3 = ("G", 1.32, INV, "n2", (8, 17, T_DIR.RIGHT))

    self.gss.run_f_debug([event1, event2, event3])
    self.gss.init_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,1),(1,2),(2,3),(4,5),(5,6),(6,7),(8,10),(10,11),
                      (12,15),(16,17),(17,19)])
  

class TestDivisions2(unittest.TestCase): # In a slightly more compex tree

  def setUp(self, genome_file=TEST_GENOME_30_6):
    params = af.prepare_genome_parameters(GENOME_PARAMS)
    #events_file = self.projbase / TEST_FOLDER2 / 'T/Events.tsv'
    events_file = TEST_FOLDER2 / 'T/Events.tsv'
    self.gss = GenomeSimulator(params, events_file, genome_file)

  def test_inversions_multiple_branches(self):

    event1 = ("G", 0.624, INV, "Root", (19, 1, T_DIR.RIGHT))
    event2 = ("G", 1.073, INV, "n1", (15, 5, T_DIR.LEFT))

    self.gss.run_f_debug([event1, event2])
    self.gss.init_divisions() 
    self.assertEqual(self.gss.initial_divisions,
                     [(0,1),(1,3),(4,6),(6,7),(8,11),(12,15),(16,19)])
  
  def test_events_for_divisions(self):

    print("***")

    #events_file = "/Users/aadavin/Desktop/Krister/Zombi/test/Events.tsv"
    
    #events = self.gss.read_genome_events_file(events_file)

    event1 = ("G", 0.024, TDUP, "Root", (2,12, T_DIR.RIGHT))
    event2 = ("G", 0.025, INV, "Root", (2, 6, T_DIR.RIGHT))
    event3 = ("G", 0.026, TDUP, "Root", (3, 7, T_DIR.RIGHT))
    event4 = ("G", 0.027, INV, "Root", (2, 11, T_DIR.RIGHT))
    event5 = ("G", 0.028, TDUP, "Root", (9, 15, T_DIR.RIGHT))
    event6 = ("G", 0.029, INV, "Root", (8, 15, T_DIR.RIGHT))

  
    #self.gss.run_f_debug([event1, event2, event3, event4, event5, event6])  

    #self.gss.run_f_debug([event1]) 
    # 
    #  

    self.gss.run_f()
    self.gss.init_divisions() 
    self.gss.redo_events_for_divisions() 
  
    for ch in self.gss.node_genomes_pieces["n1"]:
         for intergene in ch.iter_intergenes():
             if len(intergene) == 0:
               print(intergene)
             for division in intergene:
                 print(intergene, len(intergene), division, len(division))

    #self.assertEqual()

if __name__ == '__main__':
    unittest.main()
