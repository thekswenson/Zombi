"""
Simple Tree class for interacting with the extant and complete tree generated
in the T phase.
"""
from typing import Iterator
import networkx as nx

from pathlib import Path
from Bio import Phylo


class Tree:
    """
    A tree class based on a DiGraph.  It knows its root and can give you paths
    from a node to the root.

    Attibutes
    ---------
    tree: nx.DiGraph
        the tree with nodes pointing to parents
    root: str
        the root node
    """
    def __init__(self, treefile: Path):
        """
        Parameters
        ----------
        treefile : Path
            the newick file containing the tree
        """
        t = Phylo.to_networkx(Phylo.read(treefile, "newick", rooted=True)) #type: ignore
        self.tree = nx.DiGraph()
        self.root = ''
        for u, v in t.edges():
            self.root = u.name
            self.tree.add_edge(u.name, v.name)

        #Find the root:
        while self.tree.in_degree(self.root) != 0:
            self.root = next(self.tree.predecessors(self.root))

        self.tree.graph['root'] = self.root
        self._leaves = set()


    def node_set(self) -> set[str]:
        """ Get the node set of the tree. """
        return set(self.tree)


    def leaf_set(self) -> set[str]:
        """ Return the leaves of the given directed tree. """
        if self._leaves:
            return self._leaves
        
        for node in self.tree.nodes():
            if self.tree.out_degree(node) == 0:
                self._leaves.add(node)

        return self._leaves


    def path_to_root(self, node: str) -> list[str]:
        """ Return a path to the root, from the given node.  """
        path = [node]
        while node != self.root:
            node = self.parent(node)
            path.append(node)

        return path


    def parent(self, node: str) -> str:
        """ Get the parent of the given node. Return empty for the root. """
        if node == self.root:
            return ''
        else:
            return next(self.tree.predecessors(node))


    def iter_edges(self) -> Iterator[tuple[str, str]]:
        for e in self.tree.edges():
            yield e


    def __iter__(self) -> Iterator[str]:
        """ Iterate over the nodes. """
        for n in self.tree:
            yield n


def get_leaves(treefile: Path) -> set[str]:
    """
    Return the leaves of the given directed tree.

    Parameters
    ----------
    tree : nx.DiGraph
        the directed tree

    Returns
    -------
    list[str]
        the list of leaf node names
    """
    return Tree(treefile).leaf_set()
