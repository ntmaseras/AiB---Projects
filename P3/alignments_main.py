import sys
from alignment import *
from parsing import parse_matrix_and_gap_cost_in_subst_matrix,read_n_fasta
import os

def read_sequences_and_subst_matrix():    
    # Prompt the user for input
    
    list_seq = []
    
    ## read the input sequences from a fasta file
    print("Enter the path of the fasta file with the sequences: ")
    seqfile = input().strip()
    while not os.path.exists(seqfile):
        print("Error: File does not exist. Please enter a valid path:")
        seqfile = input().strip()
    list_seq = read_n_fasta(seqfile)
    
        
    ## load the substitution matrix
    print("Do you want to use the default settings for the substitution matrix and gap penalties? (y/n)")
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
            print("Invalid input. Please enter y or n.")
        
    substitution_matrix = parse_matrix_and_gap_cost_in_subst_matrix(matrix_path)
        
    return list_seq,substitution_matrix


def main():
    while True:
        
        print("Multiple alignment: do you want to align sequences with the exact method (e) or approx (a)?")
        type_alig = input().strip()
        while type_alig.lower() != 'e' and type_alig.lower() != 'a':
            print("Error. Please introduce a valid input (e/a)")
            type_alig = input().strip()
        
        
        list_of_seqs,substitution_matrix = read_sequences_and_subst_matrix()
    
        print("Do you want to see an optimal alignment and store it in a file? (y/n)")
        while True:
            see_alignment = input().strip().lower()
            if see_alignment == 'y' or see_alignment == 'n':
                break
            else: 
                print("Invalid input. Please enter Y or N.")
                
        
        ## exact solution
        if type_alig.lower() == 'e':
            score, aligned_sequences = alignment_of_3_seqs(list_of_seqs, substitution_matrix)
            if see_alignment == 'y':
                output_file = 'output/generated_alignment_exact.fasta'
                for seq in aligned_sequences:
                    print(seq)
                save_sequences_as_fasta(output_file,aligned_sequences,'')
                print("--Alignment saved in ",output_file,"--")
        ## approx solution
        elif type_alig.lower() == 'a':
            M = two_approx_algorithm_for_MSA(list_of_seqs, substitution_matrix)
            score = sp_score(M,substitution_matrix)
            if see_alignment == 'y':
                output_file = 'output/generated_alignment_approx.fasta'
                print_alignment(M)
                print_alignment(M,output_file)
                print("--Alignment saved in ",output_file,"--")
            
        print("The optimal score: ", score)
        print("Do you want to align more other sequences? (Y/N)")
        while True:
            answer = input().strip().lower()
            if answer == 'y':
                break
            elif answer == 'n':
                print('Vi ses!!')
                return
            else:
                print("Invalid input. Please enter Y or N.")
if __name__ == '__main__':
    main()