from Bio import SeqIO

def fasta_seq(input_file):
    sequence = ''
    with open(input_file,'r') as f:
        for i in SeqIO.parse(f,'fasta'): sequence = i.seq
    return str(sequence)


## generate substitution matrix and initialize gap costs 
def parse_matrix_and_gap(input_file):
    substitution_matrix = {}
    with open(input_file, 'r') as f:
        lines = f.readlines()

        # Extract keys and gap  -> special case for affine (gap extent)
        gap_cost = int(lines[0][0])
        first_line = lines[0].strip()
        
        keys = [r[0] for r in lines[1:]]
        for line in lines[1:]:
            row = line.strip().split()

            key = row[0]
            substitution_matrix[key] = {}

            for i, val in enumerate(row[1:]):
                substitution_matrix[key][keys[i]] = float(val)

    if len(first_line)> 1:
            gap_extent = int(first_line[-1])
            return substitution_matrix, [gap_cost,gap_extent]
    else:
        return substitution_matrix, gap_cost
import numpy as np

def getValues(d):    
    matrix = [[v for v in inner_dict.values()] for inner_dict in d.values()]
    keys = {k: i for i, k in enumerate(d.keys())}
    return matrix, keys


