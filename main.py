def backtrack_exact(seq1,seq2,seq3,sub_matrix,GAPCOST,alignment_matrix):
	n, m , o = len(seq1), len(seq2), len(seq3)
	v = '-'+ seq1
	w = '-' + seq2
	z = '-' + seq3
	backtrack_matrix = np.empty(n+1,m+1,o+1)
	for i in range(n):
		for j in range(m):
			for k in range(o):
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
			i, j, k -= 1
		elif backtrack_matrix[i, j, k] == 1:
			v_alig = v[i] + v_alig
			w_alig = w[j] + w_alig
			z_alig = '-' + z_alig
			i, j -= 1
		elif backtrack_matrix[i, j, k] == 2:
			v_alig = v[i] + v_alig
			w_alig = '-' + w_alig
			z_alig = z[k] + z_alig
			i,k -= 1
		elif backtrack_matrix[i, j, k] == 3:
			v_alig = '-' + v_alig
			w_alig = w[j] + w_alig
			z_alig = z[k] + z_alig	
			w,k -= 1
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


				


