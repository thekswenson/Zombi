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

      event0 = ("G", 1.21, "D", "n1", (1, 9, RIGHT))
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

      
      #self.gss.run_f_debug([event1,event2])
      #self.gss.run_f_debug([event1,  event4, event9])  
      #self.gss.obtain_divisions()       
      #self.gss.obtain_events_for_divisions()

      #for ch in self.gss.all_genomes_second["n1"]:
      #    ch.print_pieces()      
   

  def test_genetrees(self):
     
      event1 = ("G", 0.028, "P", "Root", (5, 9, RIGHT))
      
      #self.gss.run_f() 
      #self.gss.obtain_divisions()       
      #self.gss.obtain_events_for_divisions()
      
      #division_trees_folder = os.path.join('test/TestDivisions1/', "Division_trees")
      
      #self.gss.write_division_trees('test/TestDivisions1/')
      #self.gss.write_division_coordinates('test/TestDivisions1/')
      #print(self.gss.all_division_families["1"].events)


class TestTranspositions(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER1, 'T/Events.tsv')
    self.gss = GenomeSimulator(params, events_file, genome_file)



class TestDuplications(unittest.TestCase):

  def setUp(self, genome_file=TEST_GENOME):
    params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
    events_file = os.path.join(TEST_FOLDER1, 'T/Events.tsv')
    self.gss = GenomeSimulator(params, events_file, genome_file)

  def test_transpositions1(self):
     
      print("TEST")
     
      event1 = ("G", 0.11, "O", "Root", (7, RIGHT))
      event2 = ("G", 1.312, "I", "n1", (13, 15, RIGHT))
      #event3 = ("G", 1.32, "T", "n1", (17, 20, 1, RIGHT, "n2"))
      #event4 = ("G", 1.33, "T", "n1", (16, 0, 1, RIGHT, "n2"))
      #event5 = ("G", 1.34, "L", "n2", (5, 32, RIGHT))
      #event6 = ("G", 1.345, "L", "n1", (15, 3, RIGHT))
      #event7 = ("G", 1.346, "T", "n2", (1, 5, 22, RIGHT, "n1"))
      #event7 = ("G", 1.346, "I", "n2", (15, 20, RIGHT))
      

      params = af.prepare_genome_parameters(af.read_parameters(GENOME_PARAMS))
      events_file = os.path.join(TEST_FOLDER1, 'T/Events.tsv')
      

      self.gss = GenomeSimulator(params, events_file, "test/30_6.gff")
      #self.gss.run_f_debug([event1, event2, event3, event4, event5, event6, event7]) 
      #self.gss.run_f_debug([event1, event2, event3, event4, event5, event6, event7]) 
      #self.gss.run_f_debug([event1, event2, event3, event4, event5, event6, event7]) 
      self.gss.run_f_debug([event1,event2, ]) 
      #
      self.gss.obtain_divisions() 
      self.gss.obtain_events_for_divisions()
     
      for i in range(10000):

        self.gss = GenomeSimulator(params, events_file, "test/30_6.gff")
        self.gss.run_f() 
         
        events = self.gss.return_all_events()
        self.gss.obtain_divisions() 
        
        
        try:
            self.gss.obtain_events_for_divisions()

        except:
        
      
          with open("./TempEvents.txt", "w") as f:
            for event in events:
              line = ["G", str(event.time), event.etype, event.lineage]
                      
              if event.etype == "P":
                  c1, c2, c3 = event.sbpL, event.sbpR, event.sbpH
                  line.append(str((c1, c2, c3)))
              elif event.etype == "O":
                  c1 = event.sbp
              elif event.etype == "T":
                c1, c2, c3, recipient, donor = event.sbpL, event.sbpR, event.receptorsbp, event.receptorlineage, event.donorlineage
                line.append(str((c1,c2,c3,recipient, donor)))

              else:
                  c1, c2 = event.sbpL, event.sbpR
                  line.append(str((c1, c2)))

              f.write("\t".join(line)+"\n")
          print("Stopping simulation")
          break

      
      #print(self.gss.all_genomes["n1"])
      
      print("****")
      for ch in self.gss.all_genomes_second["n1"]:
          ch.print_pieces() 
      print("****")
      
      #for ch in self.gss.all_genomes["n2"]:
      #    for gene, intergene in zip(ch.genes, ch.intergenes):
      #        print("Gene", gene.total_flanking, gene.gene_family)
      #        print("Intergene", intergene.total_flanking)
      print("^^^")
      print("^^^")
      
      


if __name__ == '__main__':
    unittest.main()
