"""
Old style unittests for testing the GenomeEvents and how they map coordinates.
"""
import unittest
from pathlib import Path

from zombi.Events import Origination, Transfer, Transposition, CoordEventTwoCuts
from zombi.Events import INV, ORIG, POS, TDUP, FER
from zombi.GenomeSimulator import GenomeSimulator
from zombi.Genomes import T_DIR
import zombi.AuxiliarFunctions as af


GENOME_PARAMS = Path('tests/GenomeParametersDivisions.tsv')
TEST_FOLDER1 = Path('tests/TestDivisions1/')
TEST_GENOME = Path('tests/30_6.gff')


class TestDivisions(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME):
    params = af.prepare_genome_parameters(GENOME_PARAMS)
    events_file = TEST_FOLDER1 / 'T/Events.tsv'

    self.gss = GenomeSimulator(params, events_file, genome_file)
  
  def test_divisions(self):

      event0 = ("G", 1.21, TDUP, "n1", (1, 9, T_DIR.RIGHT))
      event1 = ("G", 1.22, INV, "n1", (15, 6, T_DIR.RIGHT))
      event2 = ("G", 1.233, INV, "n1", (11, 1, T_DIR.RIGHT))
      event3 = ("G", 1.234, INV, "n1", (8, 17, T_DIR.RIGHT))

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
      #        print("Gene", gene.total_flanking, gene.family)
      #        print("Intergene", intergene.total_flanking)

      #print("**")

  def test_divisions2(self):

      event0 = ("G", 0.028, INV, "Root", (5, 9, T_DIR.RIGHT))
      event1 = ("G", 0.21, INV, "Root", (19, 1, T_DIR.RIGHT))
      event2 = ("G", 0.37, INV, "Root", (8, 12, T_DIR.RIGHT))
      event3 = ("G", 0.40, INV, "Root", (10, 6, T_DIR.RIGHT))
      event4 = ("G", 0.44, INV, "Root", (19, 4, T_DIR.RIGHT))
      event5 = ("G", 0.48, INV, "Root", (5, 17, T_DIR.RIGHT))
      event6 = ("G", 0.53, INV, "Root", (18, 14, T_DIR.RIGHT))
      event7 = ("G", 0.644, INV, "Root", (11, 2, T_DIR.RIGHT))
      event8 = ("G", 0.74, INV, "Root", (14, 2, T_DIR.RIGHT))
      event9 = ("G", 0.748, INV, "Root", (12, 15, T_DIR.RIGHT))
      event10 = ("G", 0.919, INV, "Root", (2, 3, T_DIR.RIGHT))
      event11 = ("G", 1.27, INV, "n2", (1, 14, T_DIR.RIGHT))
      event12 = ("G", 1.42, INV, "n2", (13, 11, T_DIR.RIGHT))
      event13 = ("G", 1.48, INV, "n2", (6, 19, T_DIR.RIGHT))
      event14 = ("G", 1.55, INV, "n2", (0, 18, T_DIR.RIGHT))
      event15 = ("G", 1.47, INV, "n1", (13, 18, T_DIR.RIGHT))

      
      #self.gss.run_f_debug([event1,event2])
      #self.gss.run_f_debug([event1,  event4, event9])  
      #self.gss.obtain_divisions()       
      #self.gss.obtain_events_for_divisions()

      #for ch in self.gss.all_genomes_second["n1"]:
      #    ch.print_pieces()      
   

  def test_genetrees(self):
     
      event1 = ("G", 0.028, POS, "Root", (5, 9, T_DIR.RIGHT))
      
      #self.gss.run_f() 
      #self.gss.obtain_divisions()       
      #self.gss.obtain_events_for_divisions()
      
      #division_trees_folder = os.path.join('test/TestDivisions1/', "Division_trees")
      
      #self.gss.write_division_trees('test/TestDivisions1/')
      #self.gss.write_division_coordinates('test/TestDivisions1/')
      #print(self.gss.all_division_families["1"].events)


class TestTranspositions(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME):
    params = af.prepare_genome_parameters(GENOME_PARAMS)
    events_file = TEST_FOLDER1 / 'T/Events.tsv'
    self.gss = GenomeSimulator(params, events_file, genome_file)



class TestDuplications(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME):
    params = af.prepare_genome_parameters(GENOME_PARAMS)
    events_file = TEST_FOLDER1 / 'T/Events.tsv'
    self.gss = GenomeSimulator(params, events_file, genome_file)

  def test_transpositions1(self):
     
      print("TEST")

      event1 = ("G", 0.11, ORIG, "Root", (7, T_DIR.RIGHT))
      event2 = ("G", 1.312, INV, "n1", (13, 15, T_DIR.RIGHT))
      #event3 = ("G", 1.32, FER, "n1", (17, 20, 1, T_DIR.RIGHT, "n2"))
      #event4 = ("G", 1.33, FER, "n1", (16, 0, 1, T_DIR.RIGHT, "n2"))
      #event5 = ("G", 1.34, LOSS, "n2", (5, 32, T_DIR.RIGHT))
      #event6 = ("G", 1.345, LOSS, "n1", (15, 3, T_DIR.RIGHT))
      #event7 = ("G", 1.346, FER, "n2", (1, 5, 22, T_DIR.RIGHT, "n1"))
      #event7 = ("G", 1.346, INV, "n2", (15, 20, T_DIR.RIGHT))


      params = af.prepare_genome_parameters(GENOME_PARAMS)
      events_file = TEST_FOLDER1 / 'T/Events.tsv'
      

      self.gss = GenomeSimulator(params, events_file, TEST_GENOME)
      #self.gss.run_f_debug([event1, event2, event3, event4, event5, event6, event7]) 
      #self.gss.run_f_debug([event1, event2, event3, event4, event5, event6, event7]) 
      #self.gss.run_f_debug([event1, event2, event3, event4, event5, event6, event7]) 
      self.gss.run_f_debug([event1,event2, ]) 
      #
      self.gss.obtain_divisions() 
      self.gss.obtain_events_for_divisions()
     
      for _ in range(10000):

        self.gss = GenomeSimulator(params, events_file, TEST_GENOME)
        self.gss.run_f() 
         
        events = self.gss.return_all_events()
        self.gss.obtain_divisions() 
        
        
        try:
          self.gss.obtain_events_for_divisions()

        except:
          with open("./TempEvents.txt", "w") as f:
            for event in events:
              line = ["G", str(event.time), event.etype, event.lineage]
                      
              if event.etype == POS:
                  assert isinstance(event, Transposition)
                  c1, c2, c3 = event.sbpL, event.sbpR, event.sbpH
                  line.append(str((c1, c2, c3)))
              elif event.etype == ORIG:
                  assert isinstance(event, Origination)
                  c1 = event.sbp
              elif event.etype == FER:
                  assert isinstance(event, Transfer)
                  c1, c2, c3, recipient, donor = event.sbpL, event.sbpR, event.receptorsbp, event.receptorlineage, event.donorlineage
                  line.append(str((c1,c2,c3,recipient, donor)))

              else:
                  assert isinstance(event, CoordEventTwoCuts)
                  c1, c2 = event.sbpL, event.sbpR
                  line.append(str((c1, c2)))

              f.write("\t".join(line)+"\n")
          print("Stopping simulation")
          break

      
      #print(self.gss.all_genomes["n1"])
      
      print("****")
      for ch in self.gss.node_genomes_pieces["n1"]:
          ch.print_pieces() 
      print("****")
      
      #for ch in self.gss.all_genomes["n2"]:
      #    for gene, intergene in zip(ch.genes, ch.intergenes):
      #        print("Gene", gene.total_flanking, gene.family)
      #        print("Intergene", intergene.total_flanking)
      print("^^^")
      print("^^^")
      
      


if __name__ == '__main__':
    unittest.main()
