"""
Rules for running the Zombi simulator.

Final simulated trees will be at directories of the form:
  SIMDIR/trees/treeparams-rep0/TP1/.../T/

Final simulated genomes will be at directories of the form:
  SIMDIR/genomes/treeparams-rep0/TP1/.../genomeparams-rep0/GP1/.../G/

Final simulated sequences will be at directories of the form:
  SIMDIR/sequences/treeparams-rep0/TP1/.../genomeparams-rep0/GP1/.../sequenceparams-rep0/SP1/.../S/

- The `TP1`, `GP1`, and `SP1` are the tree, genome, and sequence parameters,
  respectively. `rep0` is the replicate number (because any set of parameters
  for each of the steps can have many replicates).
- The parameter directories (e.g. `TP1`) are named by the zombi parameters from
  the config.yaml file, and have a minus `-` separating each parameter name from
  its value.

Note that, to prevent the simulations from accidentally being rerun, there
are lock files (e.g. lock_G.flag) that are saved in the simulation directory. To
rerun the simulation, remove the lock file.
(The reason for this is that different project should share the same
simulations so that the results are comparable.)
"""
import os
import shutil

from pathlib import Path

from zombi.snakemake.parameters import modParams, PDirNames
from zombi.snakemake.parameters import getTreeParams, getGenomeParams
from zombi.snakemake.parameters import getSequenceParams, zombiSeqParamDirs
from zombi.snakemake.parameters import zombiTreeParamDirs, zombiGenomeParamDirs
from zombi.snakemake.parameters import DEFAULTTREECONFIG, DEFAULTGENOMECONFIG
from zombi.snakemake.parameters import DEFAULTSEQCONFIG

#Config File Globals:
SIMDIR = str(Path(config.get('SIMDIR', 'simulations')))  #Remove trailing slash
MAX_THREADS = config.get('MAX_THREADS', 1)
TREPS = int(config.get('TREPS', 1))
TREPS_L = list(range(TREPS))
GREPS = int(config.get('GREPS', 1))
GREPS_L = list(range(GREPS))
SREPS = int(config.get('SREPS', 1))
SREPS_L = list(range(SREPS))


#Access to config values
ZOMBI_P = config['ZOMBI']
ZOMBI_SEQP = ZOMBI_P['SEQUENCE']
ZOMBI_GENP = ZOMBI_P['GENOME']
ZOMBI_TREEP = ZOMBI_P['SPECIESTREE']

TPARAMS = PDirNames.TPARAMS.value
GPARAMS = PDirNames.GPARAMS.value
SPARAMS = PDirNames.SPARAMS.value

if not os.path.exists(DEFAULTTREECONFIG):
  raise FileNotFoundError(f'Installation problem: "{DEFAULTTREECONFIG}" not found.')
if not os.path.exists(DEFAULTGENOMECONFIG):
  raise FileNotFoundError(f'Installation problem: "{DEFAULTGENOMECONFIG}" not found.')
if not os.path.exists(DEFAULTSEQCONFIG):
  raise FileNotFoundError(f'Installation problem: "{DEFAULTSEQCONFIG}" not found.')


#### WILDCARD CONSTRAINTS ####

wildcard_constraints:
  project="[^/]+"  #Prevent `project` from matching across directories


# For the All rule
#_______________________________________________________________________________

def buildAllTargetList(wildcards):
  """
  Build the list of ultimate targets based on the settings.
  """
  files = []
  files += expand(expand(SIMDIR + '/sequences/{tparams}{gparams}{sparams}S/Genes',
                         tparams=zombiTreeParamDirs(ZOMBI_TREEP, DEFAULTTREECONFIG),
                         gparams=zombiGenomeParamDirs(ZOMBI_GENP, DEFAULTGENOMECONFIG),
                         sparams=zombiSeqParamDirs(ZOMBI_SEQP, DEFAULTSEQCONFIG)),
                  trep=TREPS_L, grep=GREPS_L, srep=SREPS_L)

  return files



# Run Zombi 
#_______________________________________________________________________________


rule zombi_run_T:
  """ Simulate trees. Remove the lock file to rerun! """
  input:
    paramfile = SIMDIR + '/trees/{tparams}/parameters/SpeciesTreeParameters.tsv',

  output:
    directory(SIMDIR + '/trees/{tparams}/T'),
    SIMDIR + '/trees/{tparams}/T/ExtantTree.nwk',

  log:
    SIMDIR + '/trees/{tparams}/logs/T.log'

  run:
    lockfile = Path(SIMDIR + f'/trees/{wildcards.tparams}/lock_T.flag')
    if lockfile.exists():
      raise WorkflowError(f'Simulation protected by lockfile. Use -t, or to '
                           f'rerun remove "{lockfile}"')

      #Run the simulation
    shell('zombi T {input.paramfile} ' + SIMDIR +
          '/trees/{wildcards.tparams} &> {log}')
    lockfile.touch()


rule zombi_run_G:
  """ Simulate genomes. Remove the lock file to rerun! """
  input:
    SIMDIR + '/genomes/{tgparams}/T',
    paramfile = SIMDIR + '/genomes/{tgparams}/parameters/GenomeParameters.tsv',

  output:
    directory(SIMDIR + '/genomes/{tgparams}/G'),
    directory(SIMDIR + '/genomes/{tgparams}/G/Genomes'),
    directory(SIMDIR + '/genomes/{tgparams}/G/Gene_families'),

  log:
    SIMDIR + '/genomes/{tgparams}/logs/G.log'

  run:
    lockfile = Path(SIMDIR + f'/genomes/{wildcards.tgparams}/lock_G.flag')
    if lockfile.exists():
      raise WorkflowError(f'Simulation protected by lockfile. Use -t, or to '
                           f'rerun remove "{lockfile}"')

      #Run the simulation
    shell('zombi G {input.paramfile} ' + SIMDIR +
          '/genomes/{wildcards.tgparams} &> {log}')
    lockfile.touch()


rule zombi_run_S:
  """ Simulate sequences. Remove the lock file to rerun! """
  input:
    SIMDIR + '/sequences/{tgsparams}/G',
    paramfile = SIMDIR + '/sequences/{tgsparams}/parameters/SequenceParameters.tsv',

  output:
    directory(SIMDIR + '/sequences/{tgsparams}/S/Genes'),

  log:
    SIMDIR + '/sequences/{tgsparams}/logs/S.log',

  threads:
    4       #Curently the S mode of Zombi only uses 4 threads maximum

  run:
    lockfile = Path(SIMDIR + f'/sequences/{wildcards.tgsparams}/lock_S.flag')
    if lockfile.exists():
      raise WorkflowError(f'Simulation protected by lockfile. Use -t, or to '
                          f'rerun remove "{lockfile}"')

      #Run the simulation
    shell('zombi S -p {threads} {input.paramfile} ' + SIMDIR +
          '/sequences/{wildcards.tgsparams} &> {log}')
    lockfile.touch()


rule zombi_link_to_T:
  """ Create a symlink to the T directory in the genomes G directory. """
  input:
    SIMDIR + '/trees/{tparams}/T',

  output:
    directory(SIMDIR + '/genomes/{tparams}/' + GPARAMS + '{gparams}/T'),

  run:
    shutil.copytree(input[0], output[0])
    #shell(f'cp -r {input[0]} {output[0]}')
    #Path(str(output)).symlink_to(f'{os.getcwd()}/{input}')
    #shell(f'ln -s {os.getcwd()}/{input[0]} {output[0]}')


rule zombi_link_to_TG:
  """ Create a symlink to the S directory in the genomes G directory. """
  input:
    treedir=SIMDIR + '/trees/{tparams}/T',
    gdir=SIMDIR + '/genomes/{tparams}/' + GPARAMS + '{gparams}G',

  output:
    treelink=directory(SIMDIR + '/sequences/{tparams}/' + GPARAMS + '{gparams}' + 
                       SPARAMS + '{sparams}/T'),
    glink=directory(SIMDIR + '/sequences/{tparams}/' + GPARAMS + '{gparams}' + 
                    SPARAMS + '{sparams}/G'),
    #genomeslink=directory(SIMDIR + '/sequences/{tparams}/' + GPARAMS +
    #                      '{gparams}' + SPARAMS + '{sparams}/G/Genomes'),

  run:
    shutil.copytree(input.treedir, output.treelink)
    try:
      shutil.copytree(input.gdir, output.glink)
    except Exception as e:
      print(f'Error copying {input.gdir} to {output.glink}: {e}')
      print(list(Path(output.glink).parent.iterdir()))
      raise
    #shell(f'cp -r {input.treedir} {output.treelink}')
    #shell(f'cp -r {input.gdir} {output.glink}')
    #Path(str(output.treelink)).symlink_to(f'{os.getcwd()}/{input.treedir}')
    #Path(str(output.glink)).symlink_to(f'{os.getcwd()}/{input.gdir}')
    #shell(f'ln -s {os.getcwd()}/{input.treedir} {output.treelink}')
    #shell(f'ln -s {os.getcwd()}/{input.gdir} {output.glink}')


# Create extra output files using zombiExporter
#____________________________________________________________________________

rule zombi_positional_orthologs:
  """
  Get the positional orthologs file from the Zombi output.
  """
  input:
    SIMDIR + '/genomes/{tgparams}/G/Gene_families',

  output:
    SIMDIR + '/genomes/{tgparams}/positional_orthologs-z_orig.json',

  shell:
    'zombiExporter po ' + SIMDIR + '/genomes/{wildcards.tgparams} {output}'


rule zombi_duplications_file:
  """
  Create the duplications file for the Zombi output.
  """
  input:
    SIMDIR + '/genomes/{tgparams}/G/Genomes',

  output:
    SIMDIR + '/genomes/{tgparams}/duplication_counts_orig.tsv',

  run:
    shell('zombiExporter dupinfo ' + SIMDIR + '/genomes/{wildcards.tgparams} {output}')


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
    SIMDIR + '/trees/{tparams}/parameters/SpeciesTreeParameters.tsv',

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
    SIMDIR + '/genomes/{tgparams}/parameters/GenomeParameters.tsv',

  run:
    outfile = Path(str(output))
    outfile.parent.mkdir(parents=True, exist_ok=True)

    #Copy the file
    shutil.copy(str(input), outfile.parent)
    #Substitute the settings
    modParams(outfile, getGenomeParams(wildcards.tgparams))


rule zombi_parameters_S:
  """
  Copy the SequenceParameters.tsv file from the resources folder and then
  substitute the setting according to the config file.
  """
  input:
    'Parameters/SequenceParameters.tsv'

  output:
    SIMDIR + '/sequences/{zparams}/parameters/SequenceParameters.tsv',

  run:
    outfile = Path(str(output))
    outfile.parent.mkdir(parents=True, exist_ok=True)

    #Copy the file
    shutil.copy(str(input), outfile.parent)
    #Substitute the settings
    modParams(outfile, getSequenceParams(wildcards.zparams))
