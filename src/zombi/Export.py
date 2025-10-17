from pathlib import Path
from collections import defaultdict as ddict
from random import choice
import re
from Bio import Phylo, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from .SequenceSimulator import node_piecesre
from .Filenames import EXTANTTREE, PRUNEDsuffix


def get_genes(genesdir: Path,
              suffix=PRUNEDsuffix) \
    -> dict[str, dict[str, dict[str, SeqRecord]]]:
    """
    Read the gene fasta records from the given directory.

    Parameters
    ----------
    genesdir : Path
        The directory containing the gene fasta files.
    suffix : str
        The suffix of the gene fasta files. Default is 'pruned.fasta'.

    Returns
    -------
    dict[str, dict[src, dict[str, SeqRecord]]]
        A mapping of gene family to gene id to genome to SeqRecord.
    """
    fam2id2genome2rec: dict[str, dict[str, dict[str, SeqRecord]]] \
        = ddict(lambda: ddict(dict))

    idre = re.compile(r'(n\d+)_(\d+)')
    files = False
    for file in genesdir.glob('*'+suffix):
        files = True
        family = file.name.replace(suffix, '')
        for record in SeqIO.parse(file, 'fasta'):
            if m := idre.match(record.id):
                fam2id2genome2rec[family][m.group(2)][m.group(1)] = record

    if not files:
        raise FileNotFoundError(f'No files found in {genesdir} with suffix '
                                f'"{suffix}".')

    return fam2id2genome2rec


def random_seq(length=100, nuc=True) -> SeqRecord:
    """
    Generate a random sequence of the given length.
    """
    if nuc:
        alphabet = 'ACGT'
    else:           #amino acids
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'

    return SeqRecord(Seq(''.join(choice(alphabet) for _ in range(length))),
                     id='random', description='')


def get_leaf_names(projdir: Path) -> list[str]:
    """
    Get the leaf names from the tree file.
    """
    treefile = projdir / 'T' / EXTANTTREE

    tree = Phylo.read(treefile, 'newick')   #type: ignore
    return [l.name for l in tree.get_terminals()]


def whole_genome_to_GFF(pieces_files: list[Path], outfile: Path,
                        leaves_only=set()):
    """
    Output the whole genomes in the `pieces_files` as a single GFF file.
    """
    pid2genome: dict[str, list[str]] = ddict(list)
    seqregions = []

        #Organize the genome pieces:
    for pieces_file in pieces_files:
        node = 'Root'
        if m:= node_piecesre.match(str(pieces_file)):
          node = m.group(1)

        if not leaves_only or node in leaves_only:
            maxcoordinate = 0
            with open(pieces_file) as f:
                f.readline()            # Header

                for line in f:
                    # ["FAMILY", "TYPE", "IDENTITY", "LENGTH", "TOTAL_LEFT", "TOTAL_RIGHT", "ORIENTATION"]
                    family, type, identity, length, tleft, tright, orientation = line.strip().split("\t")

                    maxcoordinate = max(maxcoordinate, int(tright))

                    assert type == "Gene" or type == 'Divi', 'Unexpected segment type: "{type}"'

                    pid = f'{type}_{family}'
                    pid2genome[pid].append(f'{node}\tZombi\t{type}\t{int(tleft)+1}\t'
                                           f'{tright}\t.\t{orientation}\t.\tID={pid}')

                seqregions.append(f'##sequence-region {node} 1 {maxcoordinate}')

        #Write them to the GFF file:
    with open(outfile, 'w') as f:
        f.write(f'##gff-version 3.1.26\n')
        for r in seqregions:
            f.write(f'{r}\n')
        for lines in pid2genome.values():
            for line in lines:
                f.write(f'{line}\n')
