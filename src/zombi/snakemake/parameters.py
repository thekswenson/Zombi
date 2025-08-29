"""
Functions dealing with parameters.
"""
import re
import sys

from itertools import product
from enum import StrEnum
from pathlib import Path

# The location of zombi.smk and the export.smk files:
share_zombi = Path(sys.prefix) / 'share/zombi'
rules = share_zombi / 'workflow/rules'
zombi_snakefile = rules / 'zombi.smk'
if not zombi_snakefile.exists():
  raise FileNotFoundError(f'Installation problem: "{zombi_snakefile}" not found.')
ZOMBI_SNAKEFILE = str(zombi_snakefile)
zombi_export_snakefile = rules / 'export.smk'
ZOMBI_EXPORT_SNAKEFILE = str(zombi_export_snakefile)

#parameters_dir = share_zombi / 'Parameters'
#if not parameters_dir.exists():
#  raise FileNotFoundError(f'Installation problem: "{parameters_dir}" not found.')
#TREECONFIG = str(parameters_dir / 'SpeciesTreeParameters.yaml')
#GENOMECONFIG = str(parameters_dir / 'GenomeParameters.yaml')
#SEQCONFIG = str(parameters_dir / 'SequenceParameters.yaml')

class PDirNames(StrEnum):
  TPARAMS = 'treeparams'
  GPARAMS = 'genomeparams'
  SPARAMS = 'sequenceparams'


# Generate Parameter Path Strings
#_______________________________________________________________________________


def zombiFullParamDirs(treeparams: dict[str, list], treeconfig: str,
                       genomeparams: dict[str, list], genomeconfig: str,
                       seqparams: dict[str, list], seqconfig: str) -> list[str]:
  """
  Generate the list of parameter directories.
  The directories will include replicate wildcards for the tree (trep), genome
  (grep), and sequence (srep) parameters.
  """
  treedirs = zombiTreeParamDirs(treeparams, treeconfig)
  genomedirs = zombiGenomeParamDirs(genomeparams, genomeconfig)
  seqdirs = zombiSeqParamDirs(seqparams, seqconfig)
  dirs = []
  for tdir, gdir, sdir in product(treedirs, genomedirs, seqdirs):
    dirs.append(f'{tdir}{gdir}{sdir}')

  return dirs


def expandZombiFullParamDirs(treeparams: dict[str, list], treeconfig: str,
                             genomeparams: dict[str, list], genomeconfig: str,
                             seqparams: dict[str, list], seqconfig: str,
                             treps: list[int],
                             greps: list[int],
                             sreps: list[int]) -> list[str]:
  """
  Get the full parameter directories, while expanding the replicate wildcards
  to all possible combinations.
  """
  alldirs = zombiFullParamDirs(treeparams, treeconfig,
                               genomeparams, genomeconfig,
                               seqparams, seqconfig)
  return [d.format(trep=t, grep=g, srep=s)
          for t, g, s in product(treps, greps, sreps)
          for d in alldirs]


def zombiFullParamStrs(treeparams: dict[str, list], treeconfig: str,
                       genomeparams: dict[str, list], genomeconfig: str,
                       seqparams: dict[str, list], seqconfig: str) -> list[str]:
  """
  Generate the list of parameter strings.  This is the `zombiFullParamDirs()`
  with slashes '/' replaced by underscores '_'.
  The names will include replicate wildcards for the tree (trep), genome
  (grep), and sequence (srep) parameters.
  """
  return [d.replace('/', '_').strip('_')
          for d in zombiFullParamDirs(treeparams, treeconfig,
                                      genomeparams, genomeconfig,
                                      seqparams, seqconfig)]


def zombiTreeParamDirs(treeparams: dict[str, list], defaultconfig: str) \
  -> list[str]:
  """
  Create the parameter directories for the given parameters.  Each value could
  be a single value or a list of values.
  The directories will include replicate wildcards for the tree (trep), genome
  (grep), and sequence (srep) parameters.
  """
  return [f'{PDirNames.TPARAMS}-rep{{trep}}/{d}'
          for d in zombiParamDirs(treeparams, defaultconfig)]


def zombiTreeParamStrs(treeparams: dict[str, list], defaultconfig: str) \
  -> list[str]:
  """
  Create the parameter strings for the given parameters.  This is the same as
  `zombiTreeParamDirs()`, but replaces slashes '/' with underscores '_'.
  The names will include replicate wildcards for the tree (trep), genome
  (grep), and sequence (srep) parameters.
  """
  return [d.replace('/', '_').strip('_')
          for d in zombiTreeParamDirs(treeparams, defaultconfig)]


def zombiGenomeParamDirs(genomeparams: dict[str, list], defaultconfig: str) \
  -> list[str]:
  """
  Create the parameter directories for the given parameters.  Each value could
  be a single value or a list of values.
  The directories will include replicate wildcards for the tree (trep), genome
  (grep), and sequence (srep) parameters.
  """
  return [f'{PDirNames.GPARAMS}-rep{{grep}}/{d}'
          for d in zombiParamDirs(genomeparams, defaultconfig)]

          
def zombiGenomeParamStrs(genomeparams: dict[str, list], defaultconfig: str) \
  -> list[str]:
  """
  Create the parameter strings for the given parameters.  This is the same as
  `zombiGenomeParamDirs()`, but replaces slashes '/' with underscores '_'.
  The names will include replicate wildcards for the tree (trep), genome
  (grep), and sequence (srep) parameters.
  """
  return [d.replace('/', '_').strip('_')
          for d in zombiGenomeParamDirs(genomeparams, defaultconfig)]


def zombiSeqParamDirs(seqparams: dict[str, list], defaultconfig: str) \
  -> list[str]:
  """
  Create the parameter directories for the given parameters.  Each value could
  be a single value or a list of values.
  The directories will include replicate wildcards for the tree (trep), genome
  (grep), and sequence (srep) parameters.
  """
  return [f'{PDirNames.SPARAMS}-rep{{srep}}/{d}'
          for d in zombiParamDirs(seqparams, defaultconfig)]


def zombiSeqParamStrs(seqparams: dict[str, list], defaultconfig: str) \
  -> list[str]:
  """
  Create the parameter strings for the given parameters.  This is the same as
  `zombiSeqParamDirs()`, but replaces slashes '/' with underscores '_'.
  The names will include replicate wildcards for the tree (trep), genome
  (grep), and sequence (srep) parameters.
  """
  return [d.replace('/', '_').strip('_')
          for d in zombiSeqParamDirs(seqparams, defaultconfig)]


def zombiParamDirs(params: dict[str, list], defaultconfig: str) -> list[str]:
  """
  Create the parameter part for the directory name.  Each parameter could be a
  single value or a list of values. The directory names look like:

    <PARAMETER_NAME>-<VALUE>/<PARAMETER_NAME>-<VALUE>/.../

  Notes
  -----
  - Parameters appear in the same order as in the default config file.
  """
  if not params:
    return ['']

  #Get the parameter order:
  keys = []
  for key in _getParamOrder(defaultconfig):
    if key in params:
      keys.append(key)

  pkeyset = set(params.keys())
  keyset = set(keys)
  if not pkeyset <= keyset:
    raise ValueError(f'Parameter {pkeyset - keyset} not in {defaultconfig}!')

  pathcomponents: list[list[str]] = [] #lists of strings, one for each parameter
  for key in keys:
    component = []
    if isinstance(params[key], list):
      for value in params[key]:
        component.append(f'{key}-{value}/')
    else:
      component.append(f'{key}-{params[key]}/')

    pathcomponents.append(component)

  return [''.join(combo) for combo in product(*pathcomponents)]


def _getParamOrder(defaultfile: str) -> list[str]:
  """
  Get a list of the parameters as they appear in the config file.
  """
  keys = []
  with open(defaultfile) as f:
    pattern = re.compile(r'(\S+)\s+(\S+)')
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue

      if m := re.search(pattern, line):
        param = m.group(1)
        if param in keys:
          raise ValueError(f'Duplicate parameter "{param}" found '
                           f'in {defaultfile}.')
        else:
          keys.append(param)

      else:
        raise ValueError(f'Unexpected line in "{defaultfile}":\n"{line}"')

  return keys


# Parameter File Modification
#_______________________________________________________________________________

def getTreeParams(paramspath: str) -> dict[str, str]:
  """
  Extract tree parameters from the `paramspath` directory string.
  """
  return getParamDict(paramspath, PDirNames.TPARAMS)

def getGenomeParams(paramspath: str) -> dict[str, str]:
  """
  Extract genome parameters from the `paramspath` directory string.
  """
  return getParamDict(paramspath, PDirNames.GPARAMS)

def getSequenceParams(paramspath: str) -> dict[str, str]:
  """
  Extract sequence parameters from the `paramspath` directory string.
  """
  return getParamDict(paramspath, PDirNames.SPARAMS)


def getParamDict(paramspath: str, pdirname: PDirNames) -> dict[str, str]:
  """
  Extract parameters from the section of the `paramspath` directory string
  starting with `pdirname`.
  """
  #Split the path based on the pdirname:
  if pdirname.value not in paramspath:
    raise ValueError(f'Parameter path "{paramspath}" missing delimiter '
                     f'directory "{pdirname.value}".')

  _, remainder = re.split(fr'{pdirname.value}-rep\d+', paramspath)
  if remainder and remainder[0] == '/':
    remainder = remainder[1:]

  #Remove the other parameters from the remainder:
  remainder = re.split(fr'{PDirNames.TPARAMS}-rep\d+', remainder)[0]
  remainder = re.split(fr'{PDirNames.GPARAMS}-rep\d+', remainder)[0]
  params = re.split(fr'{PDirNames.SPARAMS}-rep\d+', remainder)[0]

  #Organize the parameters by delimiter:
  if not params:
    return {}
  elif params[-1] == '/':
    params = params[:-1]

  try:
    return _splitParams(params)

  except ValueError as e:
    message = (f'We were extracting the "{pdirname.value}" parameters from '
               f'"{paramspath}".\n')
    raise ValueError(f'\nError extracting parameters from "{params}":  {e}'
                     f'\n{message}')


def _splitParams(params: str) -> dict[str, str]:
  """ Split the parameters string into a dictionary of key-value pairs. """
  paramdict = {}
  for param in params.split('/'):
    if not param:
      continue
    if '-' not in param:
      raise ValueError(f'Parameter "{param}" missing delimiter "-".')

    key, value = param.split('-')
    paramdict[key] = value

  return paramdict


def modParams(configfile: str, paramdict: dict[str, str]):
  """
  Modify the parameters in the given config file, according to the given
  parameter dictionary.
  """
  #Read the file as a single string:
  with open(configfile) as f:
    c = f.read()

  for key, value in paramdict.items():
    c = _subConfig(c, rf'(^\s*{key}\s+)\S+', value)

  with open(configfile, 'w') as f:
    f.write(c)


def _subConfig(config: str, pattern: str, value) -> str:
  """
  Modifies a pattern in the given config, if `value` is not None.
  """
  if value is None:
    return config

  m = re.search(pattern, config, flags=re.MULTILINE)
  if not m:
    raise ValueError(f'Pattern {pattern} not found in config.')

  #Verify that the pattern matches exactly once:
  if re.search(pattern, config[m.end():], flags=re.MULTILINE):
    raise ValueError(f'Multiple matches for {pattern} in parameters file.')

  return re.sub(pattern, lambda m: f'{m.group(1)}{value}', config)


def getParamValues(configfile: str) -> dict[str, str]:
  """
  Get a dictionary of default parameters from the configfile.
  """

  param2val = {}
  with open(configfile) as f:
    for line in f:
      if not line.strip() or line.startswith('#'):
        continue

      pattern = r'^(\S+)\s+(\S+)'
      m = re.search(pattern, line)
      if not m:
        raise ValueError(f'Unexpected line format in {configfile}:\n'
                         f'{line.strip()}')

      param2val[m.group(1)] = m.group(2)

  return param2val


# OLD Configfile Modification that are very specific and fragile
#_______________________________________________________________________________


## These lists define the ordering of the simulation parameters:
#TREE_PARAMS = [
#  'SPECIATION',
#  'EXTINCTION',
#  'STOPPING_RULE',
#  'TOTAL_TIME',
#  'TOTAL_LINEAGES',
#  'MIN_LINEAGES',
#  'MAX_LINEAGES',
#  'SCALE_TREE',
#  'TURNOVER',
#  'LINEAGE_PROFILE',
#  'MASS_EXTINCTION',
#  'SHIFT_SPECIATION_RATE_FREQUENCY',
#  'NUM_SPECIATION_RATE_CATEGORIES',
#  'BASE_SPECIATION',
#  'SHIFT_EXTINCTION_RATE_FREQUENCY',
#  'NUM_EXTINCTION_RATE_CATEGORIES',
#  'BASE_EXTINCTION',
#  'SEED'
#]
#GENOME_PARAMS = [
#  'TANDEMDUP',
#  'DUPLICATION',
#  'TRANSFER',
#  'LOSS',
#  'INVERSION',
#  'TRANSPOSITION',
#  'ORIGINATION',
#  'TANDEMDUP_EXTENSION',
#  'DUPLICATION_EXTENSION',
#  'TRANSFER_EXTENSION',
#  'LOSS_EXTENSION',
#  'INVERSION_EXTENSION',
#  'TRANSPOSITION_EXTENSION',
#  'REPLACEMENT_TRANSFER',
#  'ASSORTATIVE_TRANSFER',
#  'ALPHA',
#  'INITIAL_GENOME_SIZE',
#  'MIN_GENOME_SIZE',
#  'GENE_LENGTH',
#  'INTERGENE_LENGTH',
#  'PSEUDOGENIZATION',
#  'RATE_FILE',
#  'SCALE_RATES',
#  'SEED'
#]
#SEQ_PARAMS = [
#  'SCALING',
#  'SEQUENCE_SIZE',
#  'SEQUENCE',
#  'N_MODEL',
#  'ST_RATE_MULTIPLIERS',
#  'GF_RATE_MULTIPLIERS',
#  'SHIFT_SUBSTITUTION_RATE',
#  'SHIFT_CATEGORIES',
#  'BASE_RATE',
#  'SIMULATE_SEQUENCE',
#  'SCALE_GENE_TREES',
#  'SEED'
#]
#
#def modZombiTreeParams(configfile: Path,
#                       specrate=None,
#                       extrate=None,
#                       stoppingrule=None,
#                       totaltime=None,
#                       numlineages=None,
#                       minlineages=None,
#                       maxlineages=None,
#                       scaletree=None,
#                       turnover=None,
#                       lineageprofile=None,
#                       massextinction=None,
#                       shiftspecratefreq=None,
#                       numspecratecategories=None,
#                       basespeciation=None,
#                       shiftextratefreq=None,
#                       numextratecategories=None,
#                       baseextinction=None,
#                       seed=None):
#  """
#  Modifies the Zombi tree parameters in the given config.
#  Use `expandTreeParams()` to get the parameters in the correct order, and to
#  fill in the defaults.
#  """
#  #Read the file as a single string:
#  with open(configfile) as f:
#    c = f.read()
#
#  c = subConfig(c, r'(SPECIATION\s+)\w:\S+', specrate)
#  c = subConfig(c, r'(EXTINCTION\s+)\w:\S+', extrate)
#  c = subConfig(c, r'(STOPPING_RULE\s+)\d', stoppingrule)
#  c = subConfig(c, r'(TOTAL_TIME\s+)\d+', totaltime)
#  c = subConfig(c, r'(TOTAL_LINEAGES\s+)\d+', numlineages)
#  c = subConfig(c, r'(MIN_LINEAGES\s+)\d+', minlineages)
#  c = subConfig(c, r'(MAX_LINEAGES\s+)\d+', maxlineages)
#  c = subConfig(c, r'(SCALE_TREE\s+)\S+', scaletree)
#  #Tp specific:
#  c = subConfig(c, r'(TURNOVER\s+)\w:\S+', turnover)
#  c = subConfig(c, r'(LINEAGE_PROFILE\s+)\S+', lineageprofile)
#  #Tm specific:
#  c = subConfig(c, r'(MASS_EXTINCTION\s+)\S+-\S+', massextinction)
#  #Ts specific:
#  c = subConfig(c, r'(SHIFT_SPECIATION_RATE_FREQUENCY\s+)\w:\S+', shiftspecratefreq)
#  c = subConfig(c, r'(NUM_SPECIATION_RATE_CATEGORIES\s+)\d+', numspecratecategories)
#  c = subConfig(c, r'(BASE_SPECIATION\s+)\w:\S+', basespeciation)
#  c = subConfig(c, r'(SHIFT_EXTINCTION_RATE_FREQUENCY\s+)\w:\S+', shiftextratefreq)
#  c = subConfig(c, r'(NUM_EXTINCTION_RATE_CATEGORIES\s+)\d+', numextratecategories)
#  c = subConfig(c, r'(BASE_EXTINCTION\s+)\w:\S+', baseextinction)
#  c = subConfig(c, r'(SEED\s+)\d+', seed)
#
#  with open(configfile, 'w') as f:
#    f.write(c)
#
#
#def modZombiGenomeParams(configfile: Path,
#                         tduprate=None,
#                         duprate=None,
#                         transferrate=None,
#                         lossrate=None,
#                         invrate=None,
#                         transporate=None,
#                         origrate=None,
#                         tduplicationext=None,
#                         duplicationext=None,
#                         transferext=None,
#                         lossext=None,
#                         inversionext=None,
#                         transpoext=None,
#                         replacementtransfer=None,
#                         assortative_transfer=None,
#                         alpha=None,
#                         genomelen=None,
#                         mingenomelen=None,
#                         genelen=None,
#                         intergenelen=None,
#                         pseudogenization=None,
#                         ratefile=None,
#                         scalerates=None,
#                         seed=None):
#  """
#  Modifies the Zombi genome parameters in the given config.
#  Use `expandGenomeParams()` to get the parameters in the correct order, and to
#  fill in the defaults.
#  """
#  #Read the file as a single string:
#  with open(configfile) as f:
#    c = f.read()
#
#  c = subConfig(c, r'(TANDEMDUP\s+)\w:\S+', tduprate)
#  c = subConfig(c, r'(DUPLICATION\s+)\w:\S+', duprate)
#  c = subConfig(c, r'(TRANSFER\s+)\w:\S+', transferrate)
#  c = subConfig(c, r'(LOSS\s+)\w:\S+', lossrate)
#  c = subConfig(c, r'(INVERSION\s+)\w:\S+', invrate)
#  c = subConfig(c, r'(TRANSPOSITION\s+)\w:\S+', transporate)
#  c = subConfig(c, r'(ORIGINATION\s+)\w:\S+', origrate)
#  c = subConfig(c, r'(TANDEMDUP_EXTENSION\s+)\w:\S+', tduplicationext)
#  c = subConfig(c, r'(DUPLICATION_EXTENSION\s+)\w:\S+', duplicationext)
#  c = subConfig(c, r'(TRANSFER_EXTENSION\s+)\w:\S+', transferext)
#  c = subConfig(c, r'(LOSS_EXTENSION\s+)\w:\S+', lossext)
#  c = subConfig(c, r'(INVERSION_EXTENSION\s+)\w:\S+', inversionext)
#  c = subConfig(c, r'(TRANSPOSITION_EXTENSION\s+)\w:\S+', transpoext)
#  c = subConfig(c, r'(REPLACEMENT_TRANSFER\s+)\S+', replacementtransfer)
#  c = subConfig(c, r'(ASSORTATIVE_TRANSFER\s+)(True|False)', assortative_transfer)
#  c = subConfig(c, r'(ALPHA\s+)\S+', alpha)
#  c = subConfig(c, r'(INITIAL_GENOME_SIZE\s+)\d+', genomelen)
#  c = subConfig(c, r'(MIN_GENOME_SIZE\s+)\d+', mingenomelen)
#  #Gf parameters:
#  c = subConfig(c, r'(GENE_LENGTH\s+)\w:\S+', genelen)
#  c = subConfig(c, r'(INTERGENE_LENGTH\s+)\d+', intergenelen)
#  c = subConfig(c, r'(PSEUDOGENIZATION\s+)\S+', pseudogenization)
#  #Gm parameters:
#  c = subConfig(c, r'(RATE_FILE\s+)(True|False)', ratefile)
#  c = subConfig(c, r'(SCALE_RATES\s+)(True|False)', scalerates)
#  c = subConfig(c, r'(SEED\s+)\S+', seed)
#
#  with open(configfile, 'w') as f:
#    f.write(c)
#
#
#def modZombiSequenceParams(configfile: Path,
#                           scaling=None,
#                           seqlen=None,
#                           seqtype=None,
#                           nmodel=None,
#                           stratemultipliers=None,
#                           gfratemultipliers=None,
#                           shiftsubrate=None,
#                           shiftcategories=None,
#                           baserate=None,
#                           simulatesequence=None,
#                           scalegenetrees=None,
#                           seed=None):
#  """
#  Modifies the Zombi sequence parameters in the given config.
#  Use `expandSequenceParams()` to get the parameters in the correct order, and
#  to fill in the defaults.
#  """
#  #Read the file as a single string:
#  with open(configfile) as f:
#    c = f.read()
#
#  c = subConfig(c, r'(SCALING\s+)\S+', scaling)
#  c = subConfig(c, r'(SEQUENCE_SIZE\s+)\d+', seqlen)
#  c = subConfig(c, r'(SEQUENCE\s+)(nucleotide|amino-acid|codon)', seqtype)
#  c = subConfig(c, r'(N_MODEL\s+)(K2P|CUSTOM)', nmodel)
#  #Su specific:
#  c = subConfig(c, r'(ST_RATE_MULTIPLIERS\s+)\w:S+;\S+', stratemultipliers)
#  c = subConfig(c, r'(GF_RATE_MULTIPLIERS\s+)\w:S+;\S+', gfratemultipliers)
#  #Ss specific:
#  c = subConfig(c, r'(SHIFT_SUBSTITUTION_RATE\s+)\d+', shiftsubrate)
#  c = subConfig(c, r'(SHIFT_CATEGORIES\s+)\d+', shiftcategories)
#  c = subConfig(c, r'(BASE_RATE\s+)\w:\S+', baserate)
#  c = subConfig(c, r'(SIMULATE_SEQUENCE\s+)\d', simulatesequence)
#  c = subConfig(c, r'(SCALE_GENE_TREES\s+)\d', scalegenetrees)
#  c = subConfig(c, r'(SEED\s+)\d+', seed)
#
#  with open(configfile, 'w') as f:
#    f.write(c)
#
#
#def expandTreeParams(zparams: str, defaultfile: str) -> list[str]:
#  """
#  Expand the tree parameters for the given zombi parameter string, filling
#  in missing values with the defaults from the `defaultfile`.
#  """
#  return expandParams(zparams, TREE_PARAMS, defaultfile)
#
#def expandGenomeParams(zparams: str, defaultfile: str) -> list[str]:
#  """
#  Expand the genome parameters for the given zombi parameter string, filling
#  in missing values with the defaults from the `defaultfile`.
#  """
#  return expandParams(zparams, GENOME_PARAMS, defaultfile)
#
#def expandSequenceParams(zparams: str, defaultfile: str) -> list[str]:
#  """
#  Expand the sequence parameters for the given zombi parameter string, filling
#  in missing values with the defaults from the `defaultfile`.
#  """
#  return expandParams(zparams, SEQ_PARAMS, defaultfile)
#
#
#def expandParams(zparams: str, params: list[str], defaultfile: str) -> list[str]:
#  """
#  Expand the parameters. Missing parameters will be filled with the default.
#  """
#  paramdirstr = seperateParams(zparams)[params]
#  defaultparams = getDefaultParams(params, defaultfile)
#  plist = []
#  for param in params:
#    if param not in paramdirstr:
#      plist.append(defaultparams[param])
#    else:
#      _, rest = paramdirstr.split(f'{param}-')[1]
#      plist.append(rest.split('/')[0])
#
#  return plist
#
#
#def seperateParams(zparams: str) -> dict[list[str], str]:
#  """
#  Map the parameter list (e.g. TREE_PARAMS) to the parameters specified in the
#  given `zparams` directory string.
#  """
#  l2pstr = {}
#  try:
#    empty, remainder = zparams.split('treeparams/')
#    if empty != '':
#      raise ValueError(f'Unable to split on "treeparams/".')
#
#    l2pstr[TREE_PARAMS], remainder = remainder.split('genomeparams/')
#    l2pstr[GENOME_PARAMS], l2pstr[SEQ_PARAMS] = remainder.split('sequenceparams/')
#
#  except ValueError:
#    raise ValueError(f'Unexpected directory structure:\n"{zparams}"')
#
#  return l2pstr
#
#
#def getDefaultParams(params: list[str], defaultfile: str) -> dict[str, str]:
#  """
#  Get a dictionary of default tree parameters from the config.
#  """
#  with open(defaultfile) as f:
#    config = f.read()
#
#  param2val = {}
#  for param in params:
#    pattern = rf'{param}\s+(\S+)'
#    m = re.search(pattern, config)
#    if not m:
#      raise ValueError(f'Missing default parameter for {param} in '
#                       f'"{defaultfile}"!')
#
#    value = m.group(1)
#    param2val[param] = value
#
#  return param2val