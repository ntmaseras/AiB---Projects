import sys
from alignment import  *
from parsing import *
import os

def main():
    
## python3 global_alignment_command.py l/a seq1.fasta seq2.fasta subst_matrix.txt output_file.fasta
    if sys.argv[1] == 'l':
        
    
        # Initialize sequences - whether from file or from command line
        if 'fasta' in sys.argv[2]:
            seq1 = fasta_seq(sys.argv[2])
            seq2 = fasta_seq(sys.argv[3])
        else:
            seq1 = sys.argv[2]
            seq2 = sys.argv[3]
            
        # Initialize substitution matrix and gap_cost
        substitution_matrix, gap_cost = parse_matrix_and_gap(sys.argv[4])

        # Get the optimal score of the linear global alignment 
        score,_,_ = global_alignment_linear(seq1.upper(), seq2.upper(), gap_cost, substitution_matrix,sys.argv[5])
        print("The optimal score is: ",score)
        
    elif sys.argv[1] == 'a':
        if 'fasta' in sys.argv[2]:
            seq1 = fasta_seq(sys.argv[2])
            seq2 = fasta_seq(sys.argv[3])
        else:
            seq1 = sys.argv[2]
            seq2 = sys.argv[3]
            
        substitution_matrix, gaps = parse_matrix_and_gap(sys.argv[4])
        alpha = gaps[0]
        beta = gaps[1]
        

        # Get the optimal score of the linear global alignment 
        score,_,_ = global_alignment_affine(seq1.upper(), seq2.upper(), alpha,beta, substitution_matrix,sys.argv[5])
        print("The optimal score is: ",score) 
    


if __name__ == '__main__':
    main()