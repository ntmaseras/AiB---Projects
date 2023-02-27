To execute locally:

run the file src/global_alignment.py for an interactive program

for a non-interactive program, run the file src/global_alignment_command.py with the following arguments:
    for linear alignment: python3 global_alignment_command.py l seq1.fasta seq2.fasta subst_matrix.txt output_file.fasta
    for affine cost alignment: python3 global_alignment_command.py a seq1.fasta seq2.fasta subst_matrix.txt output_file.fasta


Both versions accept input sequences as strings (from the dictionary initialized by the matrix - by default A C T G) or fasta files. To parse succesfully the substitution matrix we need it in a text file in the following format:
    Linear                 Affine
5                 |     5 5
A  0  5  2  5     |     A  0  5  2  5
C  5  0  5  2     |     C  5  0  5  2  
G  2  5  0  5     |     G  2  5  0  5
T  5  2  5  0     |     T  5  2  5  0

Where the first line is the gap cost.

