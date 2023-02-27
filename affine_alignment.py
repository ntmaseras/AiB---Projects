import sys
from parsing import fasta_seq,parse_matrix_and_gap
from alignment import global_alignment_affine

def main():
    
    if 'fasta' in sys.argv[1]:
        seq1 = fasta_seq(sys.argv[1])
        seq2 = fasta_seq(sys.argv[2])
    else:
        seq1 = sys.argv[1]
        seq2 = sys.argv[2]
        
    substitution_matrix, gaps = parse_matrix_and_gap(sys.argv[3])
    alpha = gaps[0]
    beta = gaps[1]
    

    # Get the optimal score of the linear global alignment 
    score,_,_ = global_alignment_affine(seq1.upper(), seq2.upper(), alpha,beta, substitution_matrix,sys.argv[4])
    print("The optimal score is: ",score)
    


if __name__ == '__main__':
    main()