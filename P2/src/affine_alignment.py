import sys
from parsing import fasta_seq,parse_matrix_and_gap
from alignment import global_alignment_affine
import os

def affine_alignment_int():
    ## interactive function for affine alignment
    while True:
        # Prompt the user for input
        print("Enter the first sequence (or path to a fasta file):")
        arg1 = input().strip()
        print("Enter the second sequence (or path to a fasta file):")
        arg2 = input().strip()

        if 'fasta' in arg1:
            seq1 = fasta_seq(arg1)
        else:
            seq1 = arg1
        if 'fasta' in arg2:
            seq2 = fasta_seq(arg2)
        else: 
            seq2 = arg2
        
        print("Do you want to use the default settings for the substitution matrix and gap penalties? (Y/N)")
        while True:
            answer = input().strip().lower()
            if answer == 'y':
                matrix_path = 'src/input/subst_matrix_affine.txt'
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
            
        substitution_matrix, gaps = parse_matrix_and_gap(matrix_path)
        alpha = gaps[0]
        beta = gaps[1]
        
        
        print("Do you want to see an optimal alignment and store it in a file? (Y/N)")
        while True:
            answer = input().strip().lower()
            ##print(f"Using as substitution matrix: {substitution_matrix} and g(k) = {alpha} + {beta}*k")
            if answer == 'y':
                output_file = 'alignment_output.fasta'
                score,_,_ = global_alignment_affine(seq1.upper(), seq2.upper(), alpha,beta, substitution_matrix,output_file)
                print("Saved in the file ",output_file)
                break
            elif answer == 'n':
                score,_,_ = global_alignment_affine(seq1.upper(), seq2.upper(), alpha,beta, substitution_matrix)
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
    