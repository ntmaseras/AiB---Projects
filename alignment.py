from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np

def print_alignment(seq1_align,seq2_align):
    match_symb = ""
    for i in range(len(seq2_align)):
        if seq1_align[i] == seq2_align[i]:
            match_symb += "|"
        elif seq1_align[i] != seq2_align[i]:
            if (seq1_align[i] == "-" or seq2_align[i] == "-"):
                match_symb += " "
            else:
                match_symb += "*"

    print(seq1_align)
    print(match_symb)
    print(seq2_align)


def backtrack(seq1,seq2,gap_cost,M):
    alignment1, alignment2 = "", ""
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and M[i][j] == M[i-1][j] + gap_cost:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
        elif j > 0 and M[i][j] == M[i][j-1] + gap_cost:
            alignment1 = "-" + alignment1
            alignment2 = seq2[j-1] + alignment2
            j -= 1
        else:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
    return alignment1,alignment2

def global_alignment_linear(seq1, seq2, gap_cost,subst_matrix,output_file = False):
    # Initialize the matrix with zeros
    M = [[0] * (len(seq2) + 1) for i in range(len(seq1) + 1)]
    
    # Initialize the first row and column
    for i in range(1, len(seq1) + 1):
        M[i][0] = M[i-1][0] + gap_cost
    for j in range(1, len(seq2) + 1):
        M[0][j] = M[0][j-1] + gap_cost
        
    
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match_score = M[i-1][j-1] + subst_matrix[seq1[i-1]][seq2[j-1]]
            delete_score = M[i-1][j] + gap_cost
            insert_score = M[i][j-1] + gap_cost
            M[i][j] = min(match_score, delete_score, insert_score)
            
    
    score = M[-1][-1]
    
    alignment1,alignment2 = '',''
    
    if output_file:
        alignment1,alignment2 = backtrack(seq1,seq2,gap_cost,M)
        seq1_rec = SeqRecord(Seq(alignment1), id="seq1_aligned",description='')
        seq2_rec = SeqRecord(Seq(alignment2), id="seq2_aligned",description='')
        SeqIO.write([seq1_rec, seq2_rec], output_file, "fasta")
        print_alignment(alignment1,alignment2)
        
    return score,alignment1,alignment2

def backtrack_affine(seq1,seq2,alpha,beta,S,subst_matr):
    
    i = len(seq1)
    j = len(seq2)
    alignment1 = ""
    alignment2 = ""

    while i > 0 or j > 0:
        if (i > 0 and j > 0) and S[i][j] == (S[i-1][j-1] + subst_matr[seq1[i-1]][seq2[j-1]]):
            # ends in subcolumn
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
        else:
            k = 1
            while True:
                if i >= k and S[i][j] == S[i-k][j] + (alpha + k * beta):
                    # ends in deletion column
                    for n in range(k):
                        alignment1 = seq1[i-1-n] + alignment1
                        alignment2 = "-" + alignment2
                    i -= k
                    break
                elif j >= k and S[i][j] == S[i][j-k] + (alpha + k * beta):
                    # ends in insertion column
                    for n in range(k):
                        alignment1 = "-" + alignment1
                        alignment2 = seq2[j-1-n] + alignment2
                    j -= k
                    break
                else:
                    k += 1
    
    return alignment1, alignment2

def global_alignment_affine(seq1, seq2, alpha, beta, subst_matrix, output_file=False):
    n = len(seq1) + 1
    m = len(seq2) + 1

    # Initialize the matrices
    S = np.full((n, m), np.nan)
    D = np.full((n, m), np.nan)
    I = np.full((n, m), np.nan)


    # Fill in the rest of the matrices
    for i in range(n):
        for j in range(m): 
            #fill D corresponding to deletions in seq1
            d1 = d2 = np.nan # NaN values 
            if i > 0 and j >= 0:
                d1 = S[i-1,j] + (alpha + beta)
            if i > 1 and j >= 0:
                d2 = D[i-1,j] + alpha
            
            if i > 0:
                D[i,j] = np.nanmin([d1,d2]) #pick min value ignoring NaN

            #fill I The matrix with insertions (here insertions in seq1)
            i1 = i2 = np.nan
            if i >= 0 and j > 0:
                i1 = S[i,j-1] + (alpha + beta)
            if i >= 0 and j > 1:
                i2 = I[i,j-1] + alpha
            if j > 0:
                I[i][j] = np.nanmin([i1,i2])

            # fill S
            s0 = s1 = s2 = s3 = np.nan
            if i==0 and j==0:
                s0 = 0
            if i > 0 and j > 0:
                s1 = S[i-1,j-1] + subst_matrix[seq1[i-1]][seq2[j-1]] #looking in the substitution matrix
            if i > 0 and j >= 0:
                s2 = D[i][j]
            if i >=0 and j > 0:
                s3 = I[i][j]
            S[i][j] = np.nanmin([s0,s1,s2,s3])
    
    # Compute the optimal score
    score = int(S[-1][-1])
   
    alignment1,alignment2 = '',''
    
    if output_file:
        alignment1,alignment2 = backtrack_affine(seq1,seq2,alpha,beta,S,subst_matrix)
        seq1_rec = SeqRecord(Seq(alignment1), id="seq1_aligned",description='')
        seq2_rec = SeqRecord(Seq(alignment2), id="seq2_aligned",description='')
        SeqIO.write([seq1_rec, seq2_rec], output_file, "fasta")
        print_alignment(alignment1,alignment2)
        
    return score,alignment1,alignment2
    
