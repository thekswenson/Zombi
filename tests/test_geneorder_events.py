"""
Test that the gene-order events written to the files are correct.
"""
from itertools import chain
import pytest
import glob
import networkx as nx
import pandas as pd

from pathlib import Path
from Bio import Phylo
from collections import defaultdict

from zombi.Filenames import COMPLETETREE, TREEEVENTS, TREELENGTHS, EXTANTTREE
from zombi.Filenames import GENEORDEREVENTSsuffix, GENOMEsuffix
from zombi.Events import LFER_F, TDUP, DUP, LFER, AFER, LOSS, INV, POS, ORIG
from zombi.Events import LFER_B 

REPS = 10

T_PARAMS = Path('Parameters/SpeciesTreeParameters.tsv')
G_PARAMS = Path('Parameters/GenomeParameters.tsv')
#S_PARAMS = Path('Parameters/SequenceParameters.tsv')
#T_PARAMS = Path('tests/SpeciesTreeParameters.tsv') #With Seed set
#G_PARAMS = Path('tests/GenomeParameters.tsv')      #With Seed set
#S_PARAMS = Path('tests/SequenceParameters.tsv')    #With Seed set
T_SMALL_PARAMS = Path('tests/SpeciesTreeParameters_small.tsv')
Gm_PARAMS = Path('tests/GenomeParameters.tsv')

@pytest.fixture(scope='session')
def projdir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project')

@pytest.fixture
def run_T(projdir, script_runner):
  """ Test the T mode of Zombi. """
  script_runner.run(['zombi', 'T', T_PARAMS, projdir])
  return projdir / 'T'

@pytest.fixture(scope='session')
def smallprojdir(tmp_path_factory) -> Path:
  return tmp_path_factory.mktemp('zombi_project_small')

@pytest.fixture
def small_T(smallprojdir, script_runner):
  """ Test the T mode of Zombi. """
  completetree = smallprojdir / 'T' / COMPLETETREE

  while True:
    script_runner.run(['zombi', 'T', '-f', T_SMALL_PARAMS, smallprojdir])
    assert completetree.exists()

    if completetree.stat().st_size > 20:
      return smallprojdir / 'T'


def test_T(run_T):
  assert (run_T / TREEEVENTS).exists()
  assert (run_T / COMPLETETREE).exists()
  assert (run_T / EXTANTTREE).exists()
  assert (run_T / TREELENGTHS).exists()


@pytest.mark.repeat(REPS)
def test_G(projdir, script_runner, run_T):
  """ Test the G mode of Zombi. """
  assert (run_T).exists(), 'There was a problem with run_T!'

  script_runner.run(['zombi', 'G', '-a', G_PARAMS, projdir])

  outdir = projdir / 'G'
  eventsdir = outdir / 'Geneorder_events_per_branch'
  assert eventsdir.exists()
  genomesdir = outdir / 'All_genomes'
  checkEventsAgainstGenomes(run_T / COMPLETETREE, eventsdir, genomesdir)


@pytest.mark.repeat(REPS)
def test_Gf(projdir, script_runner, run_T):
  """ Test the G mode of Zombi. """
  assert (run_T).exists(), 'There was a problem with run_T!'

  script_runner.run(['zombi', 'Gf', '-a', G_PARAMS, projdir])

  outdir = projdir / 'G'
  eventsdir = outdir / 'Geneorder_events_per_branch'
  assert eventsdir.exists()
  genomesdir = outdir / 'All_genomes'
  checkEventsAgainstGenomes(run_T / COMPLETETREE, eventsdir, genomesdir)


@pytest.mark.repeat(REPS)
def test_Gm(projdir, script_runner, run_T):
  """ Test the G mode of Zombi. """
  assert (run_T).exists(), 'There was a problem with run_T!'

  script_runner.run(['zombi', 'Gm', '-a', G_PARAMS, projdir])

  outdir = projdir / 'G'
  eventsdir = outdir / 'Geneorder_events_per_branch'
  assert eventsdir.exists()
  genomesdir = outdir / 'All_genomes'
  checkEventsAgainstGenomes(run_T / COMPLETETREE, eventsdir, genomesdir)


@pytest.mark.repeat(REPS)
def test_Gu(projdir, script_runner, run_T):
  """ Test the G mode of Zombi. """
  assert (run_T).exists(), 'There was a problem with run_T!'

  script_runner.run(['zombi', 'Gu', '-a', G_PARAMS, projdir])

  outdir = projdir / 'G'
  eventsdir = outdir / 'Geneorder_events_per_branch'
  assert eventsdir.exists()
  genomesdir = outdir / 'All_genomes'
  checkEventsAgainstGenomes(run_T / COMPLETETREE, eventsdir, genomesdir)






#_______________________________________________________________________________
# Functions

def checkEventsAgainstGenomes(treefile: Path, eventsdir: Path, genomesdir: Path):
  """
  Check the the gene-order events correctly reproduce the genomes in the
  genomesdir.
  """
  tree: nx.DiGraph = Phylo.to_networkx(Phylo.read(treefile, "newick")) #type: ignore

  for parent, child in tree.edges():
    eventfile = eventsdir / f'{child}{GENEORDEREVENTSsuffix}'

    genome = getLastGenome(f'{genomesdir}/{parent}-')
    df = pd.read_csv(eventfile, sep='\t')
    #Organize events by time:
    events = defaultdict(list)
    for data in df.itertuples():
      events[data.TIME].append((data.EVENT, data.BREAKPOINTS, data.LENGTH))

    for i, (time, eventlist) in enumerate(sorted(events.items())):
      assert len(eventlist) <= 2, f'Too many events at time {time} in {eventfile}!'
      event, breakpoints, length = eventlist[0]
      if len(eventlist) > 1:
        assert event == LOSS, f'L can be the only composite event (with AT) in {eventfile}!'

      if event == TDUP:
        genome = doTandemDup(genome, breakpoints)
      elif event == DUP:
        genome = doDuplication(genome, breakpoints)
      elif event == LFER or event == LFER_F or event == LFER_B:
        pass                        #leaving transfer does not change chromosome
      elif event == AFER:           #regular (non-replacement) transfer
        genome = doArrivingTransfer(genome, int(breakpoints), time, eventsdir,
                                    genomesdir)
      elif event == LOSS:
        genome = doLoss(genome, breakpoints)
        if len(eventlist) == 2:     #replacement transfer
          assert eventlist[1][0] == AFER, f'Expected {AFER} after loss in {eventfile}!'
          genome = doArrivingTransfer(genome, int(eventlist[1][1]), time,
                                      eventsdir, genomesdir)
      elif event == INV:
        genome = doInversion(genome, breakpoints)
      elif event == POS:
        genome = doTransposition(genome, breakpoints)
      elif event == ORIG:
        genome = doOrigination(genome, int(breakpoints), length)
      else:
        raise ValueError(f'Unknown event type {event} in {eventfile}!')

      genomefile = genomesdir / f'{child}-{i}{GENOMEsuffix}'
      expected_genome = getGenome(genomefile)

      assert len(genome) == len(expected_genome), (f'Genome length mismatch '
        f'after event {event} at time {time} in {eventfile}, compared to '
        f'{genomefile}!')

      for i, (gene, egene) in enumerate(zip(genome, expected_genome)):
        if gene == 'X':
          genome[i] = egene  #Originated genes are marked as X
          gene = egene

        assert gene == egene, (f'Genome mismatch at time {time} in {eventfile},'
                               f' following event {event} at breakpoints '
                               f'{breakpoints}, compared to {genomefile}.')


def getLastGenome(fileprefix: str) -> list[str]:
  """
  Get the genome from the filename with with the highest sequence number for the
  given path prefix.
  """
  files = glob.glob(f'{fileprefix}*')
  num2file = {}
  for file in files:
    suffix = file.replace(fileprefix, '')
    num = int(suffix.replace(GENOMEsuffix, ''))
    num2file[num] = file

  return getGenome(num2file[max(num2file.keys())])


def getGenome(tsvfile: Path) -> list[str]:
  """
  Get the genome from the given filename.
  """
  #Use pandas to read the TSV file
  df = pd.read_csv(tsvfile, sep='\t')

  assert 'GENE_FAMILY' in df.columns, f'No GENE_FAMILY column in {tsvfile}!'
  assert 'ORIENTATION' in df.columns, f'No ORIENTATION column in {tsvfile}!'
  return [f'{sign}{gene}'
          for gene, sign in zip(df['GENE_FAMILY'], df['ORIENTATION'])]


def doTandemDup(genome: list[str], breakpoints: str) -> list[str]:
  """
  Perform a tandem duplication on the given genome at the given breakpoints.
  Breakpoints are given as "start,end", 0-based, inclusive.
  """
  start, end = map(int, breakpoints.split(','))
  segment = getSegment(genome, start, end)

  return genome[:end+1] + segment + genome[end+1:]


def doDuplication(genome: list[str], breakpoints: str) -> list[str]:
  """
  Perform a duplication on the given genome at the given breakpoints.
  Breakpoints are given as "start,end,here", 0-based, inclusive.
  The duplicated segment is inserted `here`.
  """
  start, end, here = map(int, breakpoints.split(','))
  segment = getSegment(genome, start, end)

  return genome[:here] + segment + genome[here:]


def doLoss(genome: list[str], breakpoints: str) -> list[str]:
  """
  Perform a loss on the given genome at the given breakpoints.
  Breakpoints are given as "start,end", 0-based, inclusive.
  """
  start, end = map(int, breakpoints.split(','))
  if start > end:
    return genome[end+1:start]

  return genome[:start] + genome[end+1:]


def doInversion(genome: list[str], breakpoints: str) -> list[str]:
  """
  Perform an inversion on the given genome at the given breakpoints.
  Breakpoints are given as "start,end", 0-based, inclusive.
  """
  start, end = map(int, breakpoints.split(','))
  if start <= end:
    positions = list(range(start, end+1))
  else:
    positions = list(chain(range(start, len(genome)), range(0, end+1)))

  revpositions = list(reversed(positions[(len(positions)+1)//2:]))
  for i, j in zip(positions[:len(positions)//2], revpositions):
    genome[i], genome[j] = genome[j], genome[i]
    genome[i] = flipSign(genome[i])
    genome[j] = flipSign(genome[j])
  
  if len(positions) % 2 == 1:
    mid = positions[len(positions)//2]
    genome[mid] = flipSign(genome[mid])
  
  return genome


def flipSigns(segment: list[str]) -> None:
  """
  Flip the signs of the genes in the given segment.
  """
  for i in range(len(segment)):
    gene = segment[i]
    segment[i] = flipSign(gene)


def flipSign(gene: str) -> str:
  """
  Flip the sign of the given gene.
  """
  if gene[0] == '+':
    return '-' + gene[1:]
  else:
    return '+' + gene[1:]


def doTransposition(genome: list[str], breakpoints: str) -> list[str]:
  """
  Perform a transposition on the given genome at the given breakpoints.
  Breakpoints are given as "start,end,here", 0-based, inclusive.
  The transposed segment is inserted `here`.
  """
  start, end, here = map(int, breakpoints.split(','))
  segment = getSegment(genome, start, end)

  #Remove the segment first
  if start <= end:
    newgenome = genome[:start] + genome[end+1:]

    #Adjust here if it is after the removed segment
    if here > start:
      here -= (end - start + 1)
  else:
    newgenome = genome[end+1:start]
    here -= end + 1

  return newgenome[:here] + segment + newgenome[here:]


def doOrigination(genome: list[str], breakpoint: int, length: int) -> list[str]:
  """
  Perform an origination on the given genome at the given breakpoints.
  Breakpoint is given as "here", 0-based, inclusive.
  The originated gene is inserted `here`.
  """
  return genome[:breakpoint] + ['X'] * length + genome[breakpoint:]


def doArrivingTransfer(genome: list[str], breakpoint: int, time: float,
                       eventsdir: Path, genomesdir: Path) -> list[str]:
  """
  Perform an arriving transfer on the given genome at the given breakpoints.
  Breakpoint is given as "here", 0-based, inclusive.  The transferred segment is
  inserted `here`, and retrieved from another file based on time.
  """
  #Find the file and line with the same time stamp
  for file in eventsdir.iterdir():
    df = pd.read_csv(file, sep='\t')
    result = df.loc[(df['TIME'] == time) & ((df['EVENT'] == LFER) |
                                            (df['EVENT'] == LFER_F) |
                                            (df['EVENT'] == LFER_B))]
    if not result.empty:
      break

  assert not result.empty, (f'No matching {LFER} event found for time {time} '
                            f'in {eventsdir}!')

  #Get the segment from the LFER event
  othergenome = getGenomeAtTime(genomesdir, file, time)
  bp1, bp2 = map(int, result['BREAKPOINTS'].item().split(','))
  segment = getSegment(othergenome, bp1, bp2)

  if result['EVENT'].item() == LFER_B:
    flipSigns(segment)
    segment.reverse()

  return genome[:breakpoint] + segment + genome[breakpoint:]


def getGenomeAtTime(genomesdir: Path, eventfile: Path, time: float) -> list[str]:
  """
  Get the genome at the given time from the given event file.
  """
  #Find the number of events before the given time:
  df = pd.read_csv(eventfile, sep='\t')
  lasttime = -1.0
  count = 0
  for row in df.itertuples():
    assert isinstance(row.TIME, float)
    if row.TIME == time:
      break

    if row.TIME > lasttime:
      count += 1
      lasttime = row.TIME

  assert row, f'No matching event found for time {time} in {eventfile}!'

  name = eventfile.name.replace(GENEORDEREVENTSsuffix, '')
  #Get the genome from the corresponding genome file
  return getGenome(genomesdir / f'{name}-{count}{GENOMEsuffix}')

  
def getSegment(genome: list[str], start: int, end: int) -> list[str]:
  """
  Get the segment from the genome between the given breakpoints.
  """
  if start > end:
    return genome[start:] + genome[:end+1]

  return genome[start:end+1]


def genomeStr(genome: list[str]) -> str:
  """
  Get a string representation of the genome for printing.
  """
  return ' '.join(genome)
