import sys
from alignment import *
from parsing import parse_matrix_and_gap_cost_in_subst_matrix
import os

def main():
    
## python3 testing_sp_exact_3.py tests/testdata_long.txt 
    
    # Open the input file for reading
    with open(sys.argv[1], 'r') as f:
        # Read the sequences and score from the file
        lines = f.readlines()
        list_of_seq = []
        for line in lines:
            if line.startswith('>'):
                # This is a sequence header line, so skip it
                continue
            else:
                # This is a sequence line, so add it to the list
                list_of_seq.append(line.strip())
        score = int(list_of_seq.pop(0))
    matrix_path = 'input/subst_matrix.txt'
    substitution_matrix = parse_matrix_and_gap_cost_in_subst_matrix(matrix_path)
    computed_score = alignment_of_3_seqs_nona(list_of_seq, substitution_matrix)
    
    if computed_score == score:
        print("Test passed")
    else:
        print("--------------Test failed------------ \nComputed score: ",computed_score,"\nGiven score: ",score,"\n-------------------------------------")
    
    
   


if __name__ == '__main__':
    main()