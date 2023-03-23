import numpy as np
import itertools as it
import pandas as pd

#Exact 3 group pairwise alignment
def alignment_of_3_seqs(list_of_seqs, subst_matrix,gap_penalty=5, print_alignment = False): #build using MSA frpm sldes (slide 19 and 20) exact
    seq1, seq2, seq3 = list_of_seqs
    t = np.zeros((len(seq1)+1, len(seq2)+1, len(seq3)+1))
    
    for i in range(len(seq1)+1):
        for j in range(len(seq2)+1):
            for k in range(len(seq3)+1):
                v1 = v2 = v3 = v4 = v5 = v6 = v7 = np.inf
                if i == 0 and j == 0 and k == 0:
                    #v0 = 0
                    t[i, j, k] = 0
                else:
                    if i > 0 and j > 0 and k > 0:
                        v1 = t[i-1, j-1, k-1] + subst_matrix[seq1[i-1]][seq2[j-1]] + subst_matrix[seq1[j-1]][seq3[k-1]] + subst_matrix[seq2[i-1]][seq3[k-1]]
                    if i > 0 and j > 0 and k >= 0:
                        v2 = t[i-1, j-1, k] + subst_matrix[seq1[i-1]][seq2[j-1]] + gap_penalty + gap_penalty
                    if i > 0 and j >= 0 and k > 0:
                        v3 = t[i-1, j, k-1] + subst_matrix[seq1[i-1]][seq3[k-1]] + gap_penalty + gap_penalty
                    if i >= 0 and j > 0 and k > 0:
                        v4 = t[i, j-1, k-1] + subst_matrix[seq2[j-1]][seq3[k-1]] + gap_penalty + gap_penalty
                    if i > 0 and j >= 0 and k >= 0:
                        v5 = t[i-1, j, k] + gap_penalty + gap_penalty
                    if i >= 0 and j > 0 and k >= 0:
                        v6 = t[i, j-1, k] + gap_penalty + gap_penalty
                    if i >= 0 and j >= 0 and k > 0:
                        v7 = t[i, j, k-1] + gap_penalty + gap_penalty
                    t[i, j, k] = min(v1, v2, v3, v4, v5, v6, v7)
    
    seq1_align, seq2_align, seq3_align = backtrack_exact(seq1,seq2,seq3,subst_matrix,5,t)
    if print_alignment:
        print(seq1_align)
        print(seq2_align)
        print(seq3_align)
    return t[len(seq1),len(seq2),len(seq3)], [seq1_align,seq2_align,seq3_align]



#MSA 2-approximation 
# first part making the pairwise global alignment
def global_alignment_linear(seq1, seq2,subst_matrix):
    # Initialize the matrix with zeros
    M = [[0] * (len(seq2) + 1) for i in range(len(seq1) + 1)]
    # Initialize the first row and column
    for i in range(1, len(seq1) + 1):
        M[i][0] = M[i-1][0] + subst_matrix['-']['-']
    for j in range(1, len(seq2) + 1):
        M[0][j] = M[0][j-1] + subst_matrix['-']['-']
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match_score = M[i-1][j-1] + subst_matrix[seq1[i-1]][seq2[j-1]]
            delete_score = M[i-1][j] + subst_matrix['-']['-']
            insert_score = M[i][j-1] + subst_matrix['-']['-']
            M[i][j] = min(match_score, delete_score, insert_score)
    return M


# Traceback and compute the alignment
def backtrack(seq1,seq2,subst_matrix,M):
    alignment1, alignment2 = "", ""
    i, j = len(seq1), len(seq2) 
    while i > 0 or j > 0:
        if i > 0 and M[i][j] == M[i-1][j] + subst_matrix['-']['-']:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
        elif j > 0 and M[i][j] == M[i][j-1] + subst_matrix['-']['-']:
            alignment1 = "-" + alignment1
            alignment2 = seq2[j-1] + alignment2
            j -= 1
        else:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
    # return list of columns of seq1 and seq2
    list_of_columns = [[i] for i in alignment1]
    for i,v in enumerate(alignment2):
        if i < len(list_of_columns):
            list_of_columns[i].append(v)
        else:
            list_of_columns[i].append("-")
            list_of_columns[i].append(v)

    return list_of_columns


def two_approx_algorithm_for_MSA(list_of_seqs, subst_matrix):
    # the first seq in list be s1 (the reference seq)
    s1 = list_of_seqs[0]
    M = None
    # make the pairwise alignment
    for i in range(1,len(list_of_seqs)):
        pairwise_matrix = global_alignment_linear(s1, list_of_seqs[i], subst_matrix)
    # make pairwise back track, corresponding to the A matrix in slides
        A = backtrack(s1, list_of_seqs[i],subst_matrix, pairwise_matrix)
    # construct matrix M
        if M is None:
            M = A
        else:
            MA = []
            i = 0
            j = 0
            while i < len(M) and j < len(A):
                # Invariant: (1) MA is a valid merge of all columns before column i in M
                # and all columns before column in A, and (2) the first row of M and A up
                # to (but not including) column i and j respectively is the same string
                # if gaps are removed.
                if M[i][0] == '-' and A[j][0] == '-':
                    # Case 1: The next column in MA is column i in M extended with the second symbol
                    # in column j in A.
                    M[i].append(A[j][1])
                    MA.append(M[i])
                    i = i + 1
                    j = j + 1
                elif M[i][0] == '-' and A[j][0] != '-':
                    # Case 2: A[j][0] is a character, so the second symbol in column j in A, A[j][1],
                    # must be in the column of MA that is the column in M where the first symbol corresponds
                    # to A[j][0]. By the invariant, this column in M is the next column in M, where the first
                    # symbol is a character, so we just moved forward in M until we find this column.
                    M[i].append('-')
                    MA.append(M[i])
                    i = i + 1
                elif M[i][0] != '-' and A[j][0] == '-':
                    # Case 3: M[i][0] is a character, so column i in M must be in the column of MA that also
                    # contains the second symbol from the column in A, where the first symbol is the character
                    # corresponding to M[i][0]. By the invariant, this column in A is the next column in A,
                    # where the first symbol is a character, so we just add columns from A to MA until we
                    # find this column.
                    c = ['-']*len(M[i])
                    c.append(A[j][1])
                    MA.append(c)
                    j = j + 1
                elif M[i][0] != '-' and A[j][0] != '-':
                    # Case 4: By the invariant the characters M[i][0] and A[j][0] are at the same position
                    # in the string spelled by the row of M and A if gaps are removed. The next column in
                    # MA is thus column i in M extended with the second symbol in column j in A.
                    M[i].append(A[j][1])
                    MA.append(M[i])
                    i = i + 1
                    j = j + 1
            if i <= len(M):
                # add the remaining coloumns of M to MA
                while i < len(M):
                    M[i].append('-')
                    MA.append(M[i]) ### THE MISTAKE WAS HERE -> MA.append(M[i].append('-')), for some reason it didn't do it well
                    i = i + 1
            if j < len(A):
                # add the remaining columns of A to MA
                k = len(MA[-1])
                while j < len(A):
                    c = ['-']*(k-1)
                    c.append(A[j][1])
                    MA.append(c)
                    j = j + 1
            M = MA
    return M


def backtrack_exact(seq1,seq2,seq3,sub_matrix,GAPCOST,alignment_matrix):
	n, m , o = len(seq1), len(seq2), len(seq3)
	v = '-'+ seq1
	w = '-' + seq2
	z = '-' + seq3
	backtrack_matrix = np.zeros((n+1,m+1,o+1)) 
	for i in range(1,n+1):
		for j in range(1,m+1):
			for k in range(1,o+1):
				all_match = alignment_matrix[i - 1,j - 1,k - 1] + sub_matrix[seq1[i - 1]][seq2[j - 1]] + sub_matrix[seq1[i - 1]][seq3[k - 1]] + sub_matrix[seq2[j - 1]][seq3[k - 1]]
				n_m_match = alignment_matrix[i - 1,j - 1,k] + sub_matrix[seq1[i - 1]][seq2[j - 1]] + GAPCOST * 2
				n_o_match = alignment_matrix[i - 1,j,k - 1] + sub_matrix[seq1[i - 1]][seq3[k - 1]] + GAPCOST * 2
				m_o_match = alignment_matrix[i,j - 1,k - 1] + sub_matrix[seq2[j - 1]][seq3[k - 1]] + GAPCOST * 2
				gap_i = alignment_matrix[i - 1,j,k] + GAPCOST * 2
				gap_j = alignment_matrix[i,j - 1,k] + GAPCOST * 2
				gap_k = alignment_matrix[i,j,k - 1] + GAPCOST * 2
				alignment_matrix[i, j, k] = min(all_match, n_m_match, n_o_match, m_o_match, gap_i, gap_j, gap_k)
				if alignment_matrix[i, j, k] == all_match:
					backtrack_matrix[i, j, k] = 0
				elif alignment_matrix[i, j, k] == n_m_match:
					backtrack_matrix[i, j, k] = 1
				elif alignment_matrix[i, j, k] == n_o_match:
					backtrack_matrix[i, j, k] = 2
				elif alignment_matrix[i, j, k] == m_o_match:
					backtrack_matrix[i, j, k] = 3
				elif alignment_matrix[i, j, k] == gap_i:
					backtrack_matrix[i, j, k] = 4
				elif alignment_matrix[i, j, k] == gap_j:
					backtrack_matrix[i, j, k] = 5
				elif alignment_matrix[i, j, k] == gap_k:
					backtrack_matrix[i, j, k] = 6
     
	i, j, k, = n, m, o
	v_alig, w_alig, z_alig = '', '', ''
	while i > 0 and j > 0 and k > 0:
		if backtrack_matrix[i,j,k]==0:
			v_alig = v[i] + v_alig
			w_alig = w[j] + w_alig
			z_alig = z[k] + z_alig
			i -= 1
			j -= 1
			k -= 1
		elif backtrack_matrix[i, j, k] == 1:
			v_alig = v[i] + v_alig
			w_alig = w[j] + w_alig
			z_alig = '-' + z_alig
			i -= 1
			j -= 1
		elif backtrack_matrix[i, j, k] == 2:
			v_alig = v[i] + v_alig
			w_alig = '-' + w_alig
			z_alig = z[k] + z_alig
			i -= 1
			k -= 1
		elif backtrack_matrix[i, j, k] == 3:
			v_alig = '-' + v_alig
			w_alig = w[j] + w_alig
			z_alig = z[k] + z_alig	
			k -= 1
			j -= 1
		elif backtrack_matrix[i,j,k] == 4:
			v_alig = v[i] + v_alig
			w_alig = '-' + w_alig
			z_alig = '-' + z_alig	
			i -= 1		
		elif backtrack_matrix[i,j,k] == 5:
			v_alig = '-' + v_alig
			w_alig = w[j] + w_alig
			z_alig = '-' + z_alig	
			j -= 1		
		elif backtrack_matrix[i,j,k] == 6:
			v_alig = '-' + v_alig
			w_alig = '-' + w_alig
			z_alig = z[k] + z_alig	
			k -= 1	
	while i > 0:
		v_alig = v[i] + v_alig
		w_alig = '-' + w_alig
		z_alig = '-' + z_alig
		i -= 1
	while j > 0:
		v_alig = '-' + v_alig
		w_alig = w[j] + w_alig
		z_alig = '-' + z_alig
		j -= 1
	while k > 0:					
		v_alig = '-' + v_alig
		w_alig = '-' + w_alig
		z_alig = z[k] + z_alig
		k -= 1

	return v_alig, w_alig, z_alig



# trying to calculate sum  of pairwise alignment score for all pairs of sequences
def unique_combinations(M: list[list]):
    """
    Precondition: `elements` does not contain duplicates.
    Postcondition: Returns unique combinations of length 2 from `elements`.

    >>> unique_combinations(["apple", "orange", "banana"])
    [("apple", "orange"), ("apple", "banana"), ("orange", "banana")]
    """
    return [list(i) for i in it.combinations([i for i,_ in enumerate(M[1])], 2)]

def sp_score(M,substitution_matrix):
    sum_score = 0
    for i in unique_combinations(M):
        for j,_ in enumerate(M):
            sum_score += substitution_matrix[M[j][i[0]]][M[j][i[1]]]
    return sum_score

    
def save_sequences_as_fasta(file_path, sequences,filename):
    with open(file_path, 'w') as fasta_file:
        for i, seq in enumerate(sequences):
            desc = filename + 'aligned'
            fasta_file.write('>seq{} {}\n'.format(i+1, desc))  # modify the description as needed
            fasta_file.write('{}\n'.format(seq))


def print_alignment(M,output_file = False,filename = ''):       
    alignment = [''.join(list(t)) for t in zip(*M)] 
    if output_file:
        save_sequences_as_fasta(output_file,alignment,filename)
    else:
        for seq in alignment:
            print(seq)


        



        



