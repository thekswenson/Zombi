"""
Tests for the command-line interface of Zombi.
We test that some of the files are created, and that the genomes created are 
consistent when using the --all-genomes flag.
"""
from zombi.Filenames import COMPLETETREE, TREEEVENTS, TREELENGTHS, EXTANTTREE
from zombi.Test import crosscheckGenomes, comparePiecesToGenomes, Filetype


def test_T(run_T):
  assert (run_T / TREEEVENTS).exists()
  assert (run_T / COMPLETETREE).exists()
  assert (run_T / EXTANTTREE).exists()
  assert (run_T / TREELENGTHS).exists()


def test_G(run_G_all):
  """ Test the G mode of Zombi. """
  assert (run_G_all / 'Genomes').exists()
  crosscheckGenomes(run_G_all)


def test_Gf(run_Gf_all):
  """ Test the Gf mode of Zombi. """
  comparePiecesToGenomes(run_Gf_all / 'Genomes')
  comparePiecesToGenomes(run_Gf_all / 'All_genomes', True)
  crosscheckGenomes(run_Gf_all)
  crosscheckGenomes(run_Gf_all, Filetype.PIECES)


def test_Gu(run_Gu):
  """ Test the Gu mode of Zombi. """
  assert (run_Gu / 'Genomes').exists()
  crosscheckGenomes(run_Gu)


def test_Gm(run_Gm):
  """ Test the Gm mode of Zombi. """
  assert (run_Gm / 'Genomes').exists()
  #crosscheckGenomes(outdir)



#_______________________________________________________________________________
# Functions
