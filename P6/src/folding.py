import random
import numpy as np


def initialize_protein_old(hp_string):
    left = {}
    right = {}
    mid = len(hp_string) // 2
    
    for pos, char in enumerate(hp_string):
        mod_pos = pos % 9
        if char == 'h':
            if mod_pos % 2 == 0 and pos < mid:
                left[mod_pos] = (pos, char)
            elif mod_pos % 2 == 1 and pos < mid:
                left[mod_pos] = (pos, char)
            elif mod_pos % 2 == 0 and pos >= mid:
                right[mod_pos] = (pos, char)
            elif mod_pos % 2 == 1 and pos >= mid:
                right[mod_pos] = (pos, char)
    
    print(left)
    print('-----------------')
    print(right)
    max_size = 0
    max_matching = ()
    
    for pos, left_tuple in left.items():
        if pos in right:
            right_tuple = right[pos]
            if left_tuple[1] == right_tuple[1]:
                matching_size = abs(left_tuple[0] - right_tuple[0]) + 1
                if matching_size > max_size:
                    max_size = matching_size
                    max_matching = (left_tuple[1], left_tuple[0], right_tuple[1], right_tuple[0])
    
    return max_matching

def initialize_protein(hp_string):
    odd = []
    even = []
    for pos, char in enumerate(hp_string):
        is_odd = pos % 2
        if char == 'h':
            if is_odd == 0:
                even.append(pos)
            else:
                odd.append(pos)
    
    print(even,odd)
    min_size = np.inf
    matching = {}
    i,j = 0,0
    for even,odd in zip(even,reversed(odd)):
        if even < odd:
            matching[even] = odd
        else:
            matching[odd] = even
        if abs(even-odd) < min_size:
            min_size = abs(even-odd)
            i,j = even,odd
    
    p = (i+j) // 2
    
    return matching,p

def score(protein):
    """
    Computes the energy of a folded protein using the HP model scoring function.
    """
    return 0


def fold_protein(d,p,n):
    
    """
    Implements the 1/4 approximation algorithm for folding a protein chain in a given
    number of iterations.
    """
    diff = 0
    res = "f"*(n-1)
    print(d)
    for key in d:
        diff = key - diff
        print(diff, key)
        if diff > 3:
            ## construct the little house
            t=(diff//2)-1
            casa = "l" + "f"*t + 2*"r" + "f"*t + "l"
            npos = key + len(casa)-1
            res= res[:key] + casa + res[npos:]
            print(res)
    print(p)
    res = res[:p] + 2 * "r" + res[p+2:]
    return res


def visualize(protein):
    """
    Generates an ASCII art representation of a folded protein chain.
    """
    return ''


def main():
    # Example usage:
    hp_string = "hphpphpphhphh"
    print(hp_string)
    matching, p = initialize_protein(hp_string.lower())
    folded = fold_protein(matching, p, len(hp_string))
    print(folded)
     
if __name__ == "__main__":
    main()