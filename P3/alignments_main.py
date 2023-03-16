import sys
from alignment import *
from parsing import fasta_seq,parse_matrix_and_gap_cost_in_subst_matrix
import os

def read_sequences_and_subst_matrix():    
    # Prompt the user for input
    
    i = 1
    list_seq = []
    
    while i<4:
        print("Enter the sequence ",i," (or path to a fasta file):")
        arg1 = input().strip()
        if 'fasta' in arg1:
            while not os.path.exists(arg1):
                print("Error: File does not exist. Please enter a valid path:")
                arg1 = input().strip()
            seq1 = fasta_seq(arg1)
            
        else:
            seq1 = arg1
        seq1 = seq1.upper()  
        list_seq.append(seq1)
        i += 1
        
    
    print("Do you want to use the default settings for the substitution matrix and gap penalties? (Y/N)")
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
            print("Invalid input. Please enter Y or N.")
        
    substitution_matrix = parse_matrix_and_gap_cost_in_subst_matrix(matrix_path)
        
    return list_seq,substitution_matrix


def main():
    
    print("Multiple alignment: do you want to align sequences with the exact method (e) or approx (a)?")
   
    type_alig = input().strip()
    while type_alig.lower() != 'e' and type_alig.lower() != 'a':
        print("Error. Please introduce a valid input (e/a)")
        type_alig = input().strip()
    
    
    list_of_seqs,substitution_matrix = read_sequences_and_subst_matrix()
    
    if type_alig.lower() == 'e':
        score = alignment_of_3_seqs(list_of_seqs, substitution_matrix)
    elif type_alig.lower() == 'a':
        M = two_approx_algorithm_for_MSA(list_of_seqs, substitution_matrix)
        score = sp_score(M,substitution_matrix)
        print(M)
    print("The score: ", score)

if __name__ == '__main__':
    main()