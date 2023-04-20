from Bio import SeqIO
from Bio import Phylo
import numpy as np

def parse_phy_file(input_file):
    matrix = []
    row_dict = {}

    with open(input_file) as f:
        lines = f.readlines()

    # extract the row labels and matrix values
    row_labels = []
    matrix_data = []
    for line in lines[1:]:
        parts = line.strip().split()
        row_labels.append(parts[0])
        matrix_data.append(list(map(float, parts[1:])))

    # create the dictionary mapping row labels to row indices
    row_dict = {label: idx for idx, label in enumerate(row_labels)}

    matrix = np.array(matrix_data)

    return matrix, row_dict



def print_tree_branch(path_tree):
    tree = Phylo.read(path_tree, 'newick')
    Phylo.draw(tree, branch_labels=lambda c: c.branch_length)
    
def print_tree(path_tree):
    tree = Phylo.read(path_tree, 'newick')
    Phylo.draw_ascii(tree)






