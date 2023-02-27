import sys
from parsing import fasta_seq,parse_matrix_and_gap
from alignment import global_alignment_linear

def main():
    ## python3 linear_alignment.py seq1.fasta seq2.fasta subst_matrix_l.txt output_file.fasta
    
    # Initialize sequences - whether from file or from command line
    if 'fasta' in sys.argv[1]:
        seq1 = fasta_seq(sys.argv[1])
        seq2 = fasta_seq(sys.argv[2])
    else:
        seq1 = sys.argv[1]
        seq2 = sys.argv[2]
        
    # Initialize substitution matrix and gap_cost
    substitution_matrix, gap_cost = parse_matrix_and_gap(sys.argv[3])

    # Get the optimal score of the linear global alignment 
    score,_,_ = global_alignment_linear(seq1.upper(), seq2.upper(), gap_cost, substitution_matrix,sys.argv[4])
    print("The optimal score is: ",score)
    


if __name__ == '__main__':
    main()