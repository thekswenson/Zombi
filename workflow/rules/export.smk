"""
These are rules built for exporting Zombi simulations to intput formats
compatible with various downstream analyses.  Currently supported are:
 - MCGP
 - FFGC

The external user of this rules library needs to ensure that OUTDIR is set in
the config file.
"""
import shutil

from pathlib import Path

from scripts.parameters import makeMCGPDataConfig
from scripts.fileconversion import exonBEDsFromBEDs

include: 'zombi.smk'

OUTDIR = config['OUTDIR']

# Process Zombi Output
#____________________________________________________________________________

rule Zombi_positional_orthologs_toproject:
  """
  Make a soft link in the inference project directory to the positional
  orthologs file.  We link to a fully qualified path so that the SIMDIR
  and OUTDIR don't need to be in the same directory.
  """
  input:
    SIMDIR + '/genomes/{zparams}/positional_orthologs-z_orig.json',

  output:
    OUTDIR + '/{project}/{zparams}/zombi/positional_orthologs-z.json',
  
  run:
    #Convert the relative path to a fully qualified path:
    inpath = Path(input[0]).resolve()
    shell("ln -s '{inpath}' '{output}'")


rule export_Zombi_to_MCGP:
  """
  Create MCGP input files in `mcgp_input/` in the simulation directory.
  """
  input:
    SIMDIR + '/sequences/{zparams}/S/Genes',
    SIMDIR + '/genomes/{zparams}/G/Genomes',
    treefile = SIMDIR + '/trees/{zparams}/T/ExtantTree.nwk',

  output:
    exonsdir = directory(SIMDIR + '/sequences/{zparams}/mcgp_input/exons'),
    dataconfig = SIMDIR + '/sequences/{zparams}/mcgp_input/input_config.yaml',
    families = SIMDIR + '/sequences/{zparams}/mcgp_input/all_gene_families.yaml',
    treefile = SIMDIR + '/sequences/{zparams}/mcgp_input/ExtantTree.nwk',

  run:
    datadir = Path(SIMDIR) / 'sequences' / wildcards.zparams / 'mcgp_input'
    #Create the BED files with families as names
    shell('zombiExporter bed ' + SIMDIR +
          '/sequences/{wildcards.zparams} {datadir}')

    #Copy the tree file to the output directory
    shutil.copy(input.treefile, output.treefile)

    #Modify the BEDs to make the exon files
    families = exonBEDsFromBEDs(datadir, Path(output.exonsdir))

    #Write the families file in JSON format
    with open(Path(output.families), 'w') as f:
      json.dump(list(families), f)

    #Make the data config file
    makeMCGPDataConfig(Path(output.dataconfig), 'ExtantTree.nwk',
                       datadir, Path(output.exonsdir))


rule export_Zombi_to_FFGC:
  """
  Create FFGC input files in `ffgc_input/` in the simulation directory.
  """
  input:
    SIMDIR + '/sequences/{zparams}/S/Genes',
    SIMDIR + '/genomes/{zparams}/G/Genomes',

  output:
    directory(SIMDIR + '/sequences/{zparams}/ffgc_input')

  shell:
    'zombiExporter ffgc ' + SIMDIR + '/sequences/{wildcards.zparams} {output}'


rule Zombi_duplicates_file_to_project:
  """
  Make a soft link in the inference project directory to the absolute path of
  the duplication counts file.
  """
  input:
    SIMDIR + '/genomes/{zparams}/duplication_counts_orig.tsv',

  output:
    OUTDIR + '/{project}/{zparams}/zombi/duplication_counts.tsv',

  run:
    #Convert the relative path to a fully qualified path:
    inpath = Path(input[0]).resolve()
    shell("ln -s '{inpath}' '{output}'")
