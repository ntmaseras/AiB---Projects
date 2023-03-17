import sys
from alignment import *
from parsing import fasta_seq,parse_matrix_and_gap_cost_in_subst_matrix
import os

def main():
    
## python3 alignments_command.py seq1.fasta seq2.fasta seq3.fasta subst_matrix.txt 
    
    # Initialize sequences - whether from file or from command line

    if 'fasta' in sys.argv[1]:
        seq1 = fasta_seq(sys.argv[1])
        seq2 = fasta_seq(sys.argv[2])
        seq3 = fasta_seq(sys.argv[3])
        
    else:
        seq1 = sys.argv[1].upper()
        seq2 = sys.argv[2].upper()
        seq3 = sys.argv[3].upper()
    
    seq1 = 'GTTCCGAAAGGCTAGCGCTAGGCGCC' #27
    seq2 = 'ATGGATTTATCTGCTCTTCG'# 21#
    seq3 = 'TGCATGCTGAAACTTCTCAACCA' # 24
    # Initialize substitution matrix and gap_cost
    substitution_matrix = parse_matrix_and_gap_cost_in_subst_matrix('input/subst_matrix.txt')
    # Get the optimal score of the linear global alignment 
    
    M = two_approx_algorithm_for_MSA([seq1,seq2,seq3], substitution_matrix)
    score = sp_score(M,substitution_matrix)
    
    print("The optimal score is: ",score)
    #print(M)
    alignment = print_alignment(M)
    


if __name__ == '__main__':
    main()