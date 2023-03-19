import sys
from alignment import *
from parsing import fasta_seq,parse_matrix_and_gap_cost_in_subst_matrix,read_n_fasta
import os

def main():

## python3 sp_exact_3.py tests/tests.fasta input/subst_matrix.txt output.fasta
    
    # Initialize sequences from file
    list_of_seqs = read_n_fasta(sys.argv[1])
    # Initialize substitution matrix and gap_cost
    substitution_matrix = parse_matrix_and_gap_cost_in_subst_matrix(sys.argv[2])
    # Get the optimal score 
    score,aligned_sequences = alignment_of_3_seqs(list_of_seqs, substitution_matrix)
    
    ## print alignment
    print('---------ALIGNMENT----------------')
    for seq in aligned_sequences:
        print(seq)
    
    ## save alignment in file
    output_file = "output/" + sys.argv[3]
    save_sequences_as_fasta(output_file,aligned_sequences,'')
    print("The optimal score is: ",score)
    print("Alignment saved in ",output_file)
    
    
   

if __name__ == '__main__':
    main()