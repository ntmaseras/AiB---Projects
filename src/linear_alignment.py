import sys
from parsing import fasta_seq,parse_matrix_and_gap
from alignment import global_alignment_linear
import os


def linear_alignment_int():
    
    while True:
        # Prompt the user for input
        print("Enter the first sequence (or path to a fasta file):")
        arg1 = input().strip()

        if 'fasta' in arg1:
            while not os.path.exists(arg1):
                print("Error: File does not exist. Please enter a valid path:")
                arg1 = input().strip()
            seq1 = fasta_seq(arg1)
        else:
            seq1 = arg1
        
        print("Enter the second sequence (or path to a fasta file):")
        arg2 = input().strip()
        if 'fasta' in arg2:
            while not os.path.exists(arg2):
                print("Error: File does not exist. Please enter a valid path:")
                arg2 = input().strip()
            seq2 = fasta_seq(arg2)
        else: 
            seq2 = arg2
        
        print("Do you want to use the default settings for the substitution matrix and gap penalties? (Y/N)")
        while True:
            answer = input().strip().lower()
            if answer == 'y':
                matrix_path = 'input/subst_matrix_linear.txt'
                break
            elif answer == 'n':
                print("Enter the path to the substitution matrix and gap penalty file:")
                matrix_path = input().strip()
                while not os.path.exists(matrix_path):
                    print("Error: File does not exist. Please enter a valid path:")
                    matrix_path = input().strip()
                break
            else:
                print("Invalid input. Please enter Y or N.")
            
        substitution_matrix, gap_cost = parse_matrix_and_gap(matrix_path)
        
        print("Do you want to see an optimal alignment and store it in a file? (Y/N)")
        while True:
            answer = input().strip().lower()
            if answer == 'y':
                output_file = 'alignment_output.fasta'
                score,_,_ = global_alignment_linear(seq1.upper(), seq2.upper(), gap_cost, substitution_matrix,output_file)
                print("Saved in the file ",output_file)
                break
            elif answer == 'n':
                score,_,_ = global_alignment_linear(seq1.upper(), seq2.upper(), gap_cost, substitution_matrix,output_file)
                break
            else: 
                print("Invalid input. Please enter Y or N.")
                
        print("The optimal score is: ",score)
        print("Do you want to align two other sequences? (Y/N)")
        while True:
            answer = input().strip().lower()
            if answer == 'y':
                break
            elif answer == 'n':
                print('Vi ses!!')
                return
            else:
                print("Invalid input. Please enter Y or N.")

'''
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
    '''