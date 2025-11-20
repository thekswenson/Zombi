"""
These are rules built for exporting Zombi simulations to intput formats
compatible with various downstream analyses.  See the other exports for
details on exporting for:
 - MCGP
 - FFGC

!!The external user of this rules library needs to ensure that OUTDIR is set in
the config file!!
"""
import shutil

from pathlib import Path

from zombi.snakemake.parameters import zombiFullParamDirs
from zombi.snakemake.parameters import expandZombiFullParamDirs
from zombi.snakemake.parameters import zombiFullParamStrs
from zombi.snakemake.parameters import expandZombiFullParamStrs
from zombi.snakemake.parameters import zombiTreeParamDirs
from zombi.snakemake.parameters import zombiTreeParamStrs
from zombi.snakemake.parameters import zombiGenomeParamDirs
from zombi.snakemake.parameters import expandZombiGenomeParamDirs
from zombi.snakemake.parameters import zombiGenomeParamStrs
from zombi.snakemake.parameters import zombiSeqParamDirs
from zombi.snakemake.parameters import zombiSeqParamStrs
from zombi.snakemake.parameters import DEFAULTTREECONFIG, DEFAULTGENOMECONFIG
from zombi.snakemake.parameters import DEFAULTSEQCONFIG, PATH_TO_RULES

include: 'zombi.smk'

OUTDIR = Path(config['OUTDIR'])


# List of application specific export files:

ZOMBI_EXPORT_MCGP_SNAKEFILE = str(PATH_TO_RULES / 'export_MCGP.smk')
ZOMBI_EXPORT_FFGC_SNAKEFILE = str(PATH_TO_RULES / 'export_FFGC.smk')


# Helpful constants for Snakefiles to import when specifying paths
# See also these in zombi.smk:
# TREPS, GREPS, SREPS, TREPS_L, GREPS_L, SREPS_L,
# ZOMBI_P, ZOMBI_TREEP, ZOMBI_GENP, ZOMBI_SEQP
#____________________________________________________________________________

# The REPS versions of these variables have replicate wildcards to be
# completed (trep, grep, srep).
ZOMBIPARAMDIRS_REPS = zombiFullParamDirs(ZOMBI_P, DEFAULTTREECONFIG,
                                         DEFAULTGENOMECONFIG, DEFAULTSEQCONFIG)
#ZOMBIPARAMDIRS included from zombi.smk
ZOMBIPARAMSTRS_REPS = zombiFullParamStrs(ZOMBI_P, DEFAULTTREECONFIG,
                                         DEFAULTGENOMECONFIG, DEFAULTSEQCONFIG)
ZOMBIPARAMSTRS = expandZombiFullParamStrs(ZOMBI_P, DEFAULTTREECONFIG,
                                          DEFAULTGENOMECONFIG, DEFAULTSEQCONFIG,
                                          TREPS_L, GREPS_L, SREPS_L)
ZOMBITREEPARAMDIRS = zombiTreeParamDirs(ZOMBI_P['TMODE'], ZOMBI_TREEP,
                                        DEFAULTTREECONFIG)
ZOMBITREEPARAMSTRS = zombiTreeParamStrs(ZOMBI_P['TMODE'], ZOMBI_TREEP,
                                        DEFAULTTREECONFIG)
ZOMBIGENOMEPARAMDIRS_REPS = zombiGenomeParamDirs(ZOMBI_P, DEFAULTTREECONFIG,
                                                 DEFAULTGENOMECONFIG)
ZOMBIGENOMEPARAMDIRS = expandZombiGenomeParamDirs(ZOMBI_P, DEFAULTTREECONFIG,
                                                  DEFAULTGENOMECONFIG,
                                                  TREPS_L, GREPS_L)
ZOMBIGENOMEPARAMSTRS = zombiGenomeParamStrs(ZOMBI_P['GMODE'], ZOMBI_GENP,
                                            DEFAULTGENOMECONFIG)
ZOMBISEQPARAMDIRS = zombiSeqParamDirs(ZOMBI_P['SMODE'], ZOMBI_SEQP,
                                      DEFAULTSEQCONFIG)
ZOMBISEQPARAMSTRS = zombiSeqParamStrs(ZOMBI_P['SMODE'], ZOMBI_SEQP,
                                      DEFAULTSEQCONFIG)


# Zombi Exports (zombiExporter output added to the export directory)
#____________________________________________________________________________

rule Zombi_export_blocks:
  """
  Use zombiExport to compute the blocks directory for each simulation.
  """
  input:
    G=SIMDIR / 'genomes/{zgproj}/G'
  output:
    directory(SIMDIR / 'genomes/{zgproj}/export/blocks')
  log:
    SIMDIR / 'genomes/{zgproj}/logs/blocks.log'
  params:
    projdir=subpath(input.G, parent=True)

  shell:
    'zombiExporter bed --no-sequences "{params.projdir}" "{output}" &> "{log}"'


rule Zombi_export_breakpoints:
  """
  Use zombiExport to compute the breaks file for each simulation.
  """
  input:
    G=SIMDIR / 'genomes/{zgproj}/G'
  output:
    breakpoints=SIMDIR / 'genomes/{zgproj}/export/breakpoints.tsv'
  log:
    SIMDIR / 'genomes/{zgproj}/logs/breakpoints.log'
  params:
    projdir=subpath(input.G, parent=True)

  shell:
    'zombiExporter breakpoints "{params.projdir}" "{output.breakpoints}" &> "{log}"'


# Process Zombi Outpu
#____________________________________________________________________________

rule Zombi_duplicates_file_to_project:
  """
  Make a soft link in the inference project directory to the absolute path of
  the duplication counts file.
  """
  input:
    SIMDIR / 'genomes/{tgparams}/duplication_counts_orig.tsv',

  output:
    OUTDIR / '{project}/{tgparams}' / SPARAMS / '{sparams}/zombi/duplication_counts.tsv',

  run:
    #Convert the relative path to a fully qualified path:
    inpath = Path(input[0]).resolve()
    shell("ln -s '{inpath}' '{output}'")


rule Zombi_positional_orthologs_toproject:
  """
  Make a soft link in the inference project directory to the positional
  orthologs file.  We link to a fully qualified path so that the SIMDIR
  and OUTDIR don't need to be in the same directory.
  """
  input:
    SIMDIR / 'genomes/{tgparams}/positional_orthologs-z_orig.json',

  output:
    OUTDIR / '{project}/{tgparams}' / SPARAMS / '{sparams}/zombi/positional_orthologs-z.json',

  run:
    #Convert the relative path to a fully qualified path:
    inpath = Path(input[0]).resolve()
    shell("ln -s '{inpath}' '{output}'")
