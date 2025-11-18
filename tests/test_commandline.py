"""
Tests for the command-line interface of Zombi.
We test that some of the files are created, and that the genomes created are 
consistent when using the --all-genomes flag.
"""
from zombi.Filenames import COMPLETETREE, TREEEVENTS, TREELENGTHS, EXTANTTREE
from zombi.Test import crosscheck_genomes, compare_pieces_to_genomes, Filetype


def test_T(run_T):
  assert (run_T / TREEEVENTS).exists()
  assert (run_T / COMPLETETREE).exists()
  assert (run_T / EXTANTTREE).exists()
  assert (run_T / TREELENGTHS).exists()


def test_G(run_G_all):
  """ Test the G mode of Zombi. """
  assert (run_G_all.G / 'Genomes').exists()
  crosscheck_genomes(run_G_all.G)


def test_Gf(run_Gf_all):
  """ Test the Gf mode of Zombi. """
  compare_pieces_to_genomes(run_Gf_all.G / 'Genomes')
  compare_pieces_to_genomes(run_Gf_all.G / 'All_genomes', True)
  crosscheck_genomes(run_Gf_all.G)
  crosscheck_genomes(run_Gf_all.G, Filetype.PIECES)


def test_Gu(run_Gu):
  """ Test the Gu mode of Zombi. """
  assert (run_Gu.G / 'Genomes').exists()
  crosscheck_genomes(run_Gu.G)


def test_Gm(run_Gm):
  """ Test the Gm mode of Zombi. """
  assert (run_Gm.G / 'Genomes').exists()
  #crosscheckGenomes(outdir)



#_______________________________________________________________________________
# Functions
