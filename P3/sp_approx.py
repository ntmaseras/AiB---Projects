import sys
from alignment import *
from parsing import fasta_seq,parse_matrix_and_gap_cost_in_subst_matrix, read_n_fasta
import os

def main():
    
## python3 sp_approx_3.py input/sequences.fasta input/subst_matrix.txt output.fasta
    
    # Initialize sequences from file
    list_of_seqs = read_n_fasta(sys.argv[1])
    # Initialize substitution matrix and gap_cost
    substitution_matrix = parse_matrix_and_gap_cost_in_subst_matrix(sys.argv[2])
    # Get the optimal score 
    M = two_approx_algorithm_for_MSA(list_of_seqs, substitution_matrix)
    score = sp_score(M,substitution_matrix)
    
    
   ## print alignment
    print('---------ALIGNMENT----------------')
    alignment = print_alignment(M)
    print("The optimal score is: ",score)
    
    output_file = "output/" + sys.argv[3]
    print_alignment(M,output_file)
    print("Alignment saved in",output_file)


if __name__ == '__main__':
    main()