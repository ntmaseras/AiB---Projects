import sys
from alignment import *
from parsing import *
import os
from msa_sp_score_3k  import compute_sp_score


def main():
    while True:
        
        print("Multiple alignment: do you want to align sequences with the exact method (e) or approx (a)?")
        type_alig = input().strip()
        while type_alig.lower() != 'e' and type_alig.lower() != 'a':
            print("Error. Please introduce a valid input (e/a)")
            type_alig = input().strip()
        
        
        list_of_seqs,substitution_matrix = read_sequences_and_subst_matrix()
        print(substitution_matrix)
    
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
            print(M)
            #score = sp_score(M,substitution_matrix)
            if see_alignment == 'y':
                output_file = 'output/generated_alignment_approx.fasta'
                print_alignment(M)
                print_alignment(M,output_file)
                print("--Alignment saved in ",output_file,"--")
                print("The score: ",compute_sp_score(output_file))
            
        #print("The optimal score: ", score)
        print("Do you want to align other sequences? (Y/N)")
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