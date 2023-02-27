import sys
from parsing import parse_matrix_and_gap
from alignment import global_alignment_linear

import time
import random
import matplotlib.pyplot as plt


def main():
    
    substitution_matrix, gap_cost = parse_matrix_and_gap(sys.argv[1])

    # Generate random sequences of varying lengths
    seq_lengths = range(100, 2000, 100)
    seqs = []
    for l in seq_lengths:
        seqs.append((''.join([random.choice(['A', 'C', 'G', 'T']) for _ in range(l)]),
                    ''.join([random.choice(['A', 'C', 'G', 'T']) for _ in range(l)])))
        
    
    times = []
    for s1, s2 in seqs:
        start_time = time.time()
        global_alignment_linear(s1, s2, gap_cost,substitution_matrix)
        end_time = time.time()
        times.append(end_time - start_time)
        print('1')

    # Plot the results
    plt.plot(seq_lengths, times)
    plt.xlabel('Sequence length')
    plt.ylabel('Running time (seconds)')
    plt.title('Running time of global_alignment_linear')
    plt.show()
    


if __name__ == '__main__':
    main()