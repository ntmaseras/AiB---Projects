import sys
from alignment import *
from parsing import *
import os
from msa_sp_score_3k  import compute_sp_score


def main():
    while True:
        
        print("--- Multiple alignment ---")
        
        
        
        list_of_seqs,substitution_matrix = read_sequences_and_subst_matrix()

                
        
        ## exact solution
        
        print("---EXACT ALIGNMENT---")
        score, aligned_sequences = alignment_of_3_seqs(list_of_seqs, substitution_matrix)
        output_file = 'output/generated_alignment_exact.fasta'
        for seq in aligned_sequences:
            print(seq)
        save_sequences_as_fasta(output_file,aligned_sequences,'')
        print("The exact score: ",compute_sp_score(output_file))
        ##print("--Alignment saved in ",output_file,"--")
        ## approx solution
        print("---APPROX ALIGNMENT---")
        M = two_approx_algorithm_for_MSA(list_of_seqs, substitution_matrix)        
        output_file = 'output/generated_alignment_approx.fasta'
        print_alignment(M)
        print_alignment(M,output_file)
        ##print("--Alignment saved in ",output_file,"--")
        print("The approx score: ",compute_sp_score(output_file))
            
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