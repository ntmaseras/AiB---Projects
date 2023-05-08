import random
import numpy as np
from hpview3k import *


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
    odds = []
    evens = []
    for pos, char in enumerate(hp_string):
        is_odd = pos % 2
        if char == 'h':
            if is_odd == 0:
                evens.append(pos)
            else:
                odds.append(pos)
    
    #print(even,odd)
    min_size = np.inf
    matching_evens = {}
    i_even,j_even = 0,0
    for even,odd in zip(evens,reversed(odds)):
        matching_evens[even] = odd
        
            
        if abs(even-odd) < min_size:
            min_size = abs(even-odd)
            i_even,j_even = even,odd
    
        
    min_size = np.inf
    i,j = 0,0
    matching_odds = {}
    for odd,even in zip(odds,reversed(evens)):

        matching_odds[odd] = even
            
        if abs(even-odd) < min_size:
            min_size = abs(even-odd)
            i,j = even,odd
    
    if len(matching_evens) > len(matching_odds):
        matching = matching_evens
        i,j = i_even, j_even
    else:
        matching = matching_odds
        
    p = (i+j) // 2
    
    return dict(sorted(matching.items())),p




def fold_protein(d,p,n):
    
    """
    Implements the 1/4 approximation algorithm for folding a protein chain in a given
    number of iterations.
    """
    
    res = "f"*(n-1)
    #print(d)
    keys = list(d.keys())
    #values = list(reversed(sorted(list(d.values()))))
    it = -1
    res = res[:p] + 2 * "R" + res[p+2:]
    prev_key = keys[0]
    print(res)
    ## build s1
    for key in d:
        
        diff = abs(key - prev_key)
        prev_key = key
        #print(key,diff,it)
        if diff > 2:
            ## construct the little house
            t = (diff//2)-1
            casa = "l" + "f"*t + 2*"r" + "f"*t + "l"
            npos = keys[it] + len(casa)-1
            res= res[:keys[it]] + casa + res[npos+1:]
        it += 1

    print(res)


    
    print('--------------------')
    print(res)
    print(len(res))

    return res


def path(hp_string):
    matching, p = initialize_protein(hp_string.lower())
    folded = fold_protein(matching, p, len(hp_string))
    
    score = try_test(hp_string,folded)
    seq = HPFold(hp_string.lower())
    path = "".join(path)
    print(path)
    seq.SetRelFold(make_relfold(path))
    seq.PrintFold()
    return score





def main():
    # Example usage:
    ## this one is OK
    hp_string = "ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh"
    print(len(hp_string))
    ## this one is NOT OK 
    #hp_string = "hphhhhhh"
    #print(hp_string)
    matching, p = initialize_protein(hp_string.lower())
    folded = fold_protein(matching, p, len(hp_string))


    
    #score = try_test(hp_string,folded)

    
    #print(folded)
     
if __name__ == "__main__":
    main()
    
    