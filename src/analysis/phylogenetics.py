from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_trees
from Bio import AlignIO, Phylo
import matplotlib.pyplot as plt
import numpy as np

class PhylogeneticAnalyzer:
    def __init__(self, alignment_path):
        """
        Initialize with a multiple sequence alignment file
        :param alignment_path: Path to alignment file (FASTA/CLUSTAL)
        """
        self.alignment = AlignIO.read(alignment_path, "fasta")
        self.tree = None
        self.distance_matrix = None

    def calculate_distances(self, model='identity'):
        """
        Calculate genetic distance matrix
        :param model: Distance model (identity, blastn, etc.)
        """
        calculator = DistanceCalculator(model)
        self.distance_matrix = calculator.get_distance(self.alignment)
        return self.distance_matrix

    def build_tree(self, method='upgma', bootstrap=100):
        """
        Construct phylogenetic tree with optional bootstrapping
        :param method: Tree construction method (upgma/nj)
        :param bootstrap: Number of bootstrap replicates
        """
        constructor = DistanceTreeConstructor()
        
        if bootstrap > 1:
            trees = bootstrap_trees(self.alignment, bootstrap, constructor, method)
            self.tree = consensus(trees)
        else:
            self.tree = constructor.nj(self.distance_matrix) if method == 'nj' \
                        else constructor.upgma(self.distance_matrix)
        
        return self.tree

    def visualize_tree(self, output_file=None):
        """Visualize phylogenetic tree with Matplotlib"""
        fig = plt.figure(figsize=(12, 8), dpi=100)
        ax = fig.add_subplot(111)
        
        Phylo.draw(self.tree, 
                  axes=ax,
                  branch_labels=lambda c: f"{c.branch_length:.3f}" if c.branch_length else "",
                  do_show=False)
        
        plt.title("Phylogenetic Tree", fontsize=14)
        if output_file:
            plt.savefig(output_file, bbox_inches='tight')
        else:
            plt.show()
        plt.close()

    def export_tree(self, output_file):
        """Export tree in Newick format"""
        Phylo.write(self.tree, output_file, "newick")

def consensus(trees):
    """Calculate majority-consensus tree from bootstrap replicates"""
    from Bio.Phylo.Consensus import _count_clades
    clade_counts = _count_clades(trees)
    tree = trees[0]
    for clade in tree.get_nonterminals():
        clade.confidence = clade_counts.get(clade, 0)/len(trees)
    return tree
