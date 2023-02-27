import sys
from parsing import parse_matrix_and_gap
from alignment import global_alignment_linear


def compare_alignments(align1, align2, ref1, ref2):
    if align1 == ref1 and align2 == ref2:
        return True
    elif align1 == ref2 and align2 == ref1:
        return True
    else:
        return False


def main():    
    ## python3 tests.py data/subst_matrix_l.txt data/tests.txt
    substitution_matrix, gap_cost = parse_matrix_and_gap(sys.argv[1])
    # Open the testing file for reading
    with open(sys.argv[2],'r') as f:
        #lines = f.readlines()
        # Read the number of tests from the first line
        num_tests = int(f.readline().strip())
        
        # Loop through each test case
        for test_num in range(num_tests):
            # Read the test type (linear or affine)
            test_type = f.readline().strip()

            # Read the two sequences to be aligned
            seq1 = f.readline().strip()
            seq2 = f.readline().strip()

            # Read the expected optimal cost of the alignment
            expected_cost = int(f.readline().strip())

            # Read the number of possible alignments
            num_alignments = int(f.readline().strip())

            # Call the global_alignment function with the input sequences and parameters
            if test_type == 'linear':
                cost,align1,align2 = global_alignment_linear(seq1.upper(), seq2.upper(), gap_cost, substitution_matrix,'out.txt')
            #elif test_type == 'affine':
            #    cost, alignments = global_alignment(seq1, seq2, affine_gap_penalty_open, affine_gap_penalty_ext, subst_matrix)

            # Check if the computed cost matches the expected cost 
            if cost != expected_cost:
                print(f"Test {test_num} failed: expected cost {expected_cost}, got {cost}")
            # Check if the optimal alignment obtained matches the one of the expected

            # Check if the computed alignments match the expected alignments
            results = []
            if num_alignments > 0:
                for alig in range(0,num_alignments):
                    ref1 = f.readline().strip()
                    ref2 = f.readline().strip()
                    results.append(compare_alignments(align1.upper(), align2.upper(), ref1.upper(), ref2.upper()))
                    
                if True not in results:
                    print(f"Test '{test_num+1}' failed: not the expected optimal alignment (check for issues tracing back)")
                
            print("Test '{test_num+1}' passed")
                

        # If all checks pass, print a success message
        print("Tests passed")



if __name__ == '__main__':
    main()