from Bio import SeqIO

def fasta_seq(input_file):
    sequence = ''
    with open(input_file,'r') as f:
        for i in SeqIO.parse(f,'fasta'): sequence = i.seq
    return str(sequence)


## generate substitution matrix and initialize gap costs 
def parse_matrix_and_gap(input_file):
    substitution_matrix = {}
    with open(input_file, 'r') as f:
        lines = f.readlines()

        # Extract keys and gap  -> special case for affine (gap extent)
        gap_cost = int(lines[0][0])
        first_line = lines[0].strip()
        
        keys = [r[0] for r in lines[1:]]
        for line in lines[1:]:
            row = line.strip().split()

            key = row[0]
            substitution_matrix[key] = {}

            for i, val in enumerate(row[1:]):
                substitution_matrix[key][keys[i]] = int(val)

    if len(first_line)> 1:
            gap_extent = int(first_line[-1])
            return substitution_matrix, [gap_cost,gap_extent]
    else:
        return substitution_matrix, gap_cost
    
def read_n_fasta(filename, n = False):
    '''
    Return the first n sequences of the fasta file. If n is not provided, return all the sequences in the file.
    '''
    sequences = [str(record.seq) for record in SeqIO.parse(filename, "fasta")]
    if n:  
        return sequences[:n]
    return sequences

def parse_matrix_and_gap_cost_in_subst_matrix(input_file):
    substitution_matrix = {}
    with open(input_file, 'r') as f:
        lines = f.readlines()

        # Extract gap cost 
        gap_cost = int(lines[0][0])

        # Extract keys and add gap symbol
        keys = [r[0] for r in lines[1:]]
        keys.append('-')

        # Parse substitution matrix
        for line in lines[1:]:
            row = line.strip().split()

            key = row[0]
            substitution_matrix[key] = {}

            for i, val in enumerate(row[1:]):
                substitution_matrix[key][keys[i]] = int(val)

            # Set gap symbol value to gap cost
            substitution_matrix[key]['-'] = gap_cost

        # Add gap symbol row
        substitution_matrix['-'] = {k: gap_cost for k in keys}

        return substitution_matrix

