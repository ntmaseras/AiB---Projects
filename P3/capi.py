import sys
from alignment import *
from parsing import *
import os
from msa_sp_score_3k  import compute_sp_score

substitution_matrix ={'A': {'A': 0, 'C': 5, 'G': 5, 'T': 5, '*': 5}, 
                    'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2, '*': 5}, 
                    'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5, '*': 5}, 
                    'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0, '*': 5},
                    '-': {'A': 5, 'C': 5, 'G': 5, 'T': 5, '*': 5},
                    '*': {'A': 5, 'C': 5, 'G': 5, 'T': 5, '*': 5}}

def backtrack(seq1,seq2, gapcost,M):
    alignment1, alignment2 = "", ""
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and M[i][j] == M[i-1][j] + gapcost:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
        elif j > 0 and M[i][j] == M[i][j-1] + gapcost:
            alignment1 = "-" + alignment1
            alignment2 = seq2[j-1] + alignment2
            j -= 1
        else:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
    return alignment1, alignment2


def pariwise_alignment(seq1, seq2,subst_matrix,gap_penalty):
    # Initialize the matrix with zeros
    M = [[0] * (len(seq2) + 1) for i in range(len(seq1) + 1)]
    # Initialize the first row and column
    for i in range(1, len(seq1) + 1):
        M[i][0] = M[i-1][0] + gap_penalty
    for j in range(1, len(seq2) + 1):
        M[0][j] = M[0][j-1] + gap_penalty
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match_score = M[i-1][j-1] + subst_matrix[seq1[i-1]][seq2[j-1]]
            delete_score = M[i-1][j] + gap_penalty
            insert_score = M[i][j-1] + gap_penalty
            M[i][j] = min(match_score, delete_score, insert_score)
    return M

def reference_alignment(seq1, seq2):
    ref_alignment = ''
    for a, b in zip(seq1, seq2):
        if a == '-' and b != '-':
            ref_alignment += b
        elif a != '-' and b == '-':
            ref_alignment += a
        elif a == b:
            ref_alignment += a
        else:
            ref_alignment += '*'
    # Handle any remaining nucleotides in the longer sequence
    if len(seq1) > len(seq2):
        ref_alignment += seq1[len(seq2):].replace('-', '')
    elif len(seq2) > len(seq1):
        ref_alignment += seq2[len(seq1):].replace('-', '')
    
    return ref_alignment

def add_row(matrix, new_row):
    # Determine the length of the new row
    new_row_length = len(new_row)
    # Determine the length of the existing rows
    existing_row_lengths = [len(row) for row in matrix]
    max_row_length = max(existing_row_lengths)
    # If the new row is shorter than the longest existing row, add dashes to the end
    if new_row_length < max_row_length:
        new_row += '-' * (max_row_length - new_row_length)
    # If the new row is longer than the longest existing row, add dashes to the end of all existing rows
    if new_row_length > max_row_length:
        for i in range(len(matrix)):
            matrix[i] += '-' * (new_row_length - max_row_length)
    matrix.append(new_row)
    return matrix

def capibara_slow_alignement(list_of_seqs,substitution_matrix,gapcost):
    #initialization
    M = pariwise_alignment(list_of_seqs[0],list_of_seqs[1],substitution_matrix,gapcost)
    v_aligned, w_aligned = backtrack(list_of_seqs[0], list_of_seqs[1], gapcost, M)
    ref_alignment = reference_alignment(v_aligned, w_aligned)
    matrix = [[i for i in v_aligned],[i for i in w_aligned]]
    print(matrix)
    #iterate over the list of sequences from 2 to the end
    for i in range(2,len(list_of_seqs)):
        # compute the pairwise alignment matrix between the reference alignment and the current sequence
        print("Aligning pairwise: ",i)
        temp = pariwise_alignment(list_of_seqs[i],ref_alignment,substitution_matrix,gapcost)
        # compute the alignment between the current sequence and the reference alignment
        v_aligned, w_aligned = backtrack(list_of_seqs[i], ref_alignment, gapcost, temp)
        # update the reference alignment
        ref_alignment = reference_alignment(v_aligned, w_aligned)
        # update the matrix
        new_row = [i for i in v_aligned]
        add_row(matrix,new_row)
    return matrix


list_of_seq = read_n_fasta('experiments/brca1-testseqs.fasta',6)
for i in [3]:
        M = capibara_slow_alignement(list_of_seq[:i],substitution_matrix,5)
        alignment_approx_outputfile = 'presentation/capi/aligned_' + str(i)
        approx_alignment = print_alignment(M,alignment_approx_outputfile) ## print alignment to file

alignment_outputfile = 'presentation/full_aligned' 
aligned_seq = [''.join(i) for i in M] 
save_sequences_as_fasta(alignment_outputfile,aligned_seq,'')