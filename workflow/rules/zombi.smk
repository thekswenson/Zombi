"""
Rules for running the Zombi simulator.

Note that, to prevent the simulations from accidentally being rerun, there
are lock files (e.g. lock_G.flag) that are saved in the simulation directory. To
rerun the simulation, remove the lock file.
(The reason for this is that different project should share the same
simulations so that the results are comparable.)
"""
import os
import shutil

from pathlib import Path

from scripts.parameters import modParams, PDirNames
from scripts.parameters import getTreeParams, getGenomeParams, getSequenceParams

#configfile: 'workflow/config/default.yaml'

#Config File Globals:
SIMDIR = str(Path(config.get('SIMDIR', 'simulations')))  #Remove trailing slash
MAX_THREADS = config.get('MAX_THREADS', 1)
REPS = int(config.get('REPS', 1))

#Access to config values
ZOMBI_P = config['ZOMBI']
ZOMBI_SEQP = ZOMBI_P['SEQUENCE']
ZOMBI_GENP = ZOMBI_P['GENOME']
ZOMBI_TREEP = ZOMBI_P['SPECIESTREE']

TPARAMS = PDirNames.TPARAMS.value
GPARAMS = PDirNames.GPARAMS.value
SPARAMS = PDirNames.SPARAMS.value


# Run Zombi 
#_______________________________________________________________________________

rule zombi_run_T:
  """ Simulate trees. Remove the lock file to rerun! """
  input:
    paramfile = SIMDIR + '/trees/{tparams}/rep{rep}/parameters/SpeciesTreeParameters.tsv',

  output:
    directory(SIMDIR + '/trees/{tparams}/rep{rep}/T'),
    SIMDIR + '/trees/{tparams}/rep{rep}/T/ExtantTree.nwk',

  log:
    SIMDIR + '/trees/{tparams}/rep{rep}/logs/T.log'

  run:
    lockfile = Path(SIMDIR + f'/trees/{wildcards.tparams}/rep{wildcards.rep}/lock_T.flag')
    if lockfile.exists():
      raise WorkflowError(f'Simulation protected by lockfile. Use -t, or to '
                           f'rerun remove "{lockfile}"')

      #Run the simulation
    shell('Zombi T {input.paramfile} ' + SIMDIR +
          '/trees/{wildcards.tparams}/rep{wildcards.rep} &> {log}')
    lockfile.touch()


rule zombi_run_G:
  """ Simulate genomes. Remove the lock file to rerun! """
  input:
    SIMDIR + '/genomes/{tgparams}/rep{rep}/T',
    paramfile = SIMDIR + '/genomes/{tgparams}/rep{rep}/parameters/GenomeParameters.tsv',

  output:
    directory(SIMDIR + '/genomes/{tgparams}/rep{rep}/G'),
    directory(SIMDIR + '/genomes/{tgparams}/rep{rep}/G/Genomes'),
    directory(SIMDIR + '/genomes/{tgparams}/rep{rep}/G/Gene_families'),

  log:
    SIMDIR + '/genomes/{tgparams}/rep{rep}/logs/G.log'

  run:
    lockfile = Path(SIMDIR + f'/genomes/{wildcards.tgparams}/rep{wildcards.rep}/lock_G.flag')
    if lockfile.exists():
      raise WorkflowError(f'Simulation protected by lockfile. Use -t, or to '
                           f'rerun remove "{lockfile}"')

      #Run the simulation
    shell('Zombi G {input.paramfile} ' + SIMDIR +
          '/genomes/{wildcards.tgparams}/rep{wildcards.rep} &> {log}')
    lockfile.touch()


rule zombi_run_S:
  """ Simulate sequences. Remove the lock file to rerun! """
  input:
    SIMDIR + '/sequences/{tgsparams}/rep{rep}/G',
    paramfile = SIMDIR + '/sequences/{tgsparams}/rep{rep}/parameters/SequenceParameters.tsv',

  output:
    directory(SIMDIR + '/sequences/{tgsparams}/rep{rep}/S/Genes'),

  log:
    SIMDIR + '/sequences/{tgsparams}/rep{rep}/logs/S.log',

  threads:
    4       #Curently the S mode of Zombi only uses 4 threads maximum

  run:
    lockfile = Path(SIMDIR + f'/sequences/{wildcards.tgsparams}/rep{wildcards.rep}/lock_S.flag')
    if lockfile.exists():
      raise WorkflowError(f'Simulation protected by lockfile. Use -t, or to '
                          f'rerun remove "{lockfile}"')

      #Run the simulation
    shell('Zombi S -p {threads} {input.paramfile} ' + SIMDIR +
          '/sequences/{wildcards.tgsparams}/rep{wildcards.rep} &> {log}')
    lockfile.touch()


rule zombi_link_to_T:
  """ Create a symlink to the T directory in the genomes G directory. """
  input:
    SIMDIR + '/trees/{tparams}/rep{rep}/T',

  output:
    directory(SIMDIR + '/genomes/{tparams}/' + GPARAMS + '{gparams}rep{rep}/T'),

  run:
    #Path(str(output)).symlink_to(f'{os.getcwd()}/{input}')
    #shell(f'ln -s {os.getcwd()}/{input[0]} {output[0]}')
    shell(f'cp -r {input[0]} {output[0]}')

rule zombi_link_to_TG:
  """ Create a symlink to the S directory in the genomes G directory. """
  input:
    treedir=SIMDIR + '/trees/{tparams}/rep{rep}/T',
    genomedir=SIMDIR + '/genomes/{tparams}/' + GPARAMS + '{gparams}rep{rep}/G',

  output:
    treelink=directory(SIMDIR + '/sequences/{tparams}/' + GPARAMS + '{gparams}' + 
                       SPARAMS + '{sparams}rep{rep}/T'),
    genomelink=directory(SIMDIR + '/sequences/{tparams}/' + GPARAMS + '{gparams}' + 
                         SPARAMS + '{sparams}rep{rep}/G'),

  run:
    #Path(str(output.treelink)).symlink_to(f'{os.getcwd()}/{input.treedir}')
    #Path(str(output.genomelink)).symlink_to(f'{os.getcwd()}/{input.genomedir}')
    #shell(f'ln -s {os.getcwd()}/{input.treedir} {output.treelink}')
    #shell(f'ln -s {os.getcwd()}/{input.genomedir} {output.genomelink}')
    shell(f'cp -r {input.treedir} {output.treelink}')
    shell(f'cp -r {input.genomedir} {output.genomelink}')


# Create extra output files using ZombiExporter
#____________________________________________________________________________

rule zombi_positional_orthologs:
  """
  Get the positional orthologs file from the Zombi output.
  """
  input:
    SIMDIR + '/genomes/{gparams}/rep{rep}/G/Gene_families',

  output:
    SIMDIR + '/genomes/{gparams}/rep{rep}/positional_orthologs-z_orig.json',

  shell:
    'ZombiExporter po ' + SIMDIR + '/{wildcards.gparams}/rep{wildcards.rep} {output}'


rule zombi_duplications_file:
  """
  Create the duplications file for the Zombi output.
  """
  input:
    SIMDIR + '/genomes/{gparams}/rep{rep}/G/Genomes',

  output:
    SIMDIR + '/genomes/{gparams}/rep{rep}/duplication_counts_orig.tsv',

  run:
    shell('ZombiExporter dupinfo ' + SIMDIR + '/{wildcards.gparams}/rep{wildcards.rep} {output}')


#rule Zombi_duplicates_file_toproject:
#  input:
#    SIMDIR + '/genomes/{zparams}/{rep}/duplication_counts_orig.tsv',
#
#  output:
#    OUTDIR + '/genomes/{project}/{zparams}/{rep}/duplication_counts.tsv',
#  
#  shell:
#    "ln -s '../../../../../{input}' '{output}'"


# Zombi Input
#____________________________________________________________________________

rule zombi_parameters_T:
  """
  Copy the SpeciesTreeParameters.tsv file from the resources folder and then
  substitute the setting according to the config file.
  """
  input:
    'Parameters/SpeciesTreeParameters.tsv'

  output:
    SIMDIR + '/trees/{tparams}/rep{rep}/parameters/SpeciesTreeParameters.tsv',

  run:
    outfile = Path(str(output))
    outfile.parent.mkdir(parents=True, exist_ok=True)

    #Copy the default settings file
    shutil.copy(str(input), outfile.parent)
    #Substitute the settings
    modParams(outfile, getTreeParams(wildcards.tparams))


rule zombi_parameters_G:
  """
  Copy the GenomeParmeters.tsv file from the resources folder and then
  substitute the setting according to the config file.
  """
  input:
    'Parameters/GenomeParameters.tsv'

  output:
    SIMDIR + '/genomes/{gparams}/rep{rep}/parameters/GenomeParameters.tsv',

  run:
    outfile = Path(str(output))
    outfile.parent.mkdir(parents=True, exist_ok=True)

    #Copy the file
    shutil.copy(str(input), outfile.parent)
    #Substitute the settings
    modParams(outfile, getGenomeParams(wildcards.gparams))


rule zombi_parameters_S:
  """
  Copy the SequenceParameters.tsv file from the resources folder and then
  substitute the setting according to the config file.
  """
  input:
    'Parameters/SequenceParameters.tsv'

  output:
    SIMDIR + '/sequences/{sparams}/rep{rep}/parameters/SequenceParameters.tsv',

  run:
    outfile = Path(str(output))
    outfile.parent.mkdir(parents=True, exist_ok=True)

    #Copy the file
    shutil.copy(str(input), outfile.parent)
    #Substitute the settings
    modParams(outfile, getSequenceParams(wildcards.sparams))
