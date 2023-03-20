from Bio import SeqIO
import sys
from alignment import *
from parsing import *
import os

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
                substitution_matrix[key][keys[i]] = int(val)

    if len(first_line)> 1:
            gap_extent = int(first_line[-1])
            return substitution_matrix, [gap_cost,gap_extent]
    else:
        return substitution_matrix, gap_cost
    
def read_n_fasta(filename, n = False):
    '''
    Return the first n sequences of the fasta file. If n is not provided, return all the sequences in the file.
    '''
    sequences = [str(record.seq) for record in SeqIO.parse(filename, "fasta")]
    if n:  
        return sequences[:n]
    return sequences

def parse_matrix_and_gap_cost_in_subst_matrix(input_file):
    substitution_matrix = {}
    with open(input_file, 'r') as f:
        lines = f.readlines()

        # Extract gap cost 
        gap_cost = int(lines[0][0])

        # Extract keys and add gap symbol
        keys = [r[0] for r in lines[1:]]
        keys.append('-')

        # Parse substitution matrix
        for line in lines[1:]:
            row = line.strip().split()

            key = row[0]
            substitution_matrix[key] = {}

            for i, val in enumerate(row[1:]):
                substitution_matrix[key][keys[i]] = int(val)

            # Set gap symbol value to gap cost
            substitution_matrix[key]['-'] = gap_cost

        # Add gap symbol row
        substitution_matrix['-'] = {k: gap_cost for k in keys}

        return substitution_matrix

def read_sequences_and_subst_matrix():    
    # Prompt the user for input
    
    list_seq = []
    
    ## read the input sequences from a fasta file
    print("Enter the path of the fasta file with the sequences: ")
    seqfile = input().strip()
    while not os.path.exists(seqfile):
        print("Error: File does not exist. Please enter a valid path:")
        seqfile = input().strip()
    list_seq = read_n_fasta(seqfile)
    
        
    ## load the substitution matrix
    print("Do you want to use the default settings for the substitution matrix and gap penalties? (y/n)")
    while True:
        answer = input().strip().lower()
        if answer == 'y':
            matrix_path = 'input/subst_matrix.txt'
            break
        elif answer == 'n':
            print("Enter the path to the substitution matrix and gap penalty file:")
            matrix_path = input().strip()
            while not os.path.exists(matrix_path):
                print("Error: File does not exist. Please enter a valid path:")
                matrix_path = input().strip()
            break
        else:
            print("Invalid input. Please enter y or n.")
        
    substitution_matrix = parse_matrix_and_gap_cost_in_subst_matrix(matrix_path)
        
    return list_seq,substitution_matrix
