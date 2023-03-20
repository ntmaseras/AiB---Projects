P3 - Bjarke, Paula, Genona

To execute locally:

run the file alignments_main.py for an interactive program

for a non-interactive program, run the file:

    - for exact solution:

    python3 sp_exact_3.py input/sequences.fasta input/subst_matrix.txt output_filename.fasta

    - for aprox: 

    python3 sp_approx.py input/sequences.fasta input/subst_matrix.txt output_filename.fasta


To parse succesfully the substitution matrix we need it in a text file in the following format:
                         

5                            

A  0  5  2  5  

C  5  0  5  2 

G  2  5  0  5 

T  5  2  5  0


Where the first line is the gap cost.

