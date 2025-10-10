"""
Export Zombi simulations to MCGP input formats.
"""
import shutil

from pathlib import Path
from pybedtools import BedTool, Interval


# Convert Zombi output to MCGP format
#____________________________________________________________________________


def exonBEDsFromBEDs(beddir: Path, outdir: Path) -> set[str]:
  """
  Convert each BED file in the bed dir to an exon BED file.

  Returns
  -------
  set[str]
    A set of the (family) names in the input BED files.
  """
  #Ensure that the output directory exists:
  outdir.mkdir(parents=True, exist_ok=True)

  familynames = set()
  for bedfile in beddir.glob('*.bed'):
    #Read the BED file:
    newfeatures = []
    for feature in BedTool(bedfile):
      familynames.add(feature.name)
      newfeatures.append(Interval(feature.chrom, feature.start, feature.end,
                                  feature.fields[6] + '-RA-E1', feature.score,
                                  feature.strand))

    #Write the new BED file:
    #Get the name of the file without the extension:
    BedTool(newfeatures).saveas(str(outdir / f'{bedfile.stem}.exon.bed'))

  return familynames


def makeMCGPDataConfig(outfile: Path, phylofile: Path, datadir: Path,
                       exonsdir: Path) -> None:
  """
  Makes a config file for the MCGP data.
  """
  with open(outfile, 'w') as f:
    f.write(f'phylogeny: {phylofile}\n')
    f.write(f'genomes:\n')
    for file in datadir.glob('*.fa'):
      f.write(f'  {file.stem}: "{file.name}"\n')

    f.write(f'homolog_families:\n')
    for file in datadir.glob('*.bed'):
      f.write(f'  {file.stem}: "{file.name}"\n')

    f.write(f'exons:\n')
    for file in exonsdir.glob('*.bed'):
      f.write(f'  {file.stem.replace(".exon", "")}: "exons/{file.name}"\n')


rule export_Zombi_to_MCGP:
  """
  Create MCGP input files in `mcgp_input/` in the simulation directory.
  """
  input:
    SIMDIR + '/sequences/{zparams}/S/Genes',
    SIMDIR + '/sequences/{zparams}/G', #/Genomes',
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
