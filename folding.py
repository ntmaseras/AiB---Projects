import random
import numpy as np
from hpview3k import *
import pandas as pd


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
    

    ## get the evens vs odds
    min_size = np.inf
    matching_evens = {}
    i_even,j_even = 0,0
    for even,odd in zip(evens,reversed(odds)):
        matching_evens[even] = odd
            
        if abs(even-odd) < min_size:
            min_size = abs(even-odd)
            i_even,j_even = even,odd
    
    
    ## get the odds vs evens
    min_size = np.inf
    i,j = 0,0
    matching_odds = {}
    for odd,even in zip(odds,reversed(evens)):

        matching_odds[odd] = even
            
        if abs(even-odd) < min_size:
            min_size = abs(even-odd)
            i,j = even,odd
    
    ## get the max 
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
    it = -1
    res = res[:p] + 2 * "R" + res[p+2:]
    prev_key = keys[0]
    ## build s1
    for key in d:
        diff = abs(key - prev_key)
        prev_key = key
        #print(key,diff,it)
        if diff > 2:
            ## construct the little house
            t = (diff//2)-1
            casa = "l" + "f"*t + 2*"r" + "f"*t + "l"
            npos = keys[it] + len(casa)
            res= res[:keys[it]] + casa + res[npos:]
        it += 1
        


    return res



def main():

    filename = 'P6/src/sequences.txt'
    sequences = {}
    with open(filename, "r") as f:
        for i, line in enumerate(f.readlines()):
            seq = line.split(' ')[0]
            score = int(line.split(' ')[1].replace('\n',''))
            sequences[i] = [seq,score]
    results = pd.DataFrame(columns=['seq','opt_score','c_score'])
    for id, [seq,score] in sequences.items():
        
        matching, p = initialize_protein(seq.lower())
        fold = fold_protein(matching, p, len(seq)) 
        
        results.loc[id,'seq'] = seq 
        seq = HPFold(seq)
        absfold = make_absfold(fold)
        relfold = make_relfold(fold)
        if len(absfold) == len(seq) - 1:
            seq.SetAbsFold(absfold)
        elif len(relfold) == len(seq) - 1:
            seq.SetRelFold(relfold)
            results.loc[id,'opt_score'] = score 
            results.loc[id,'c_score'] = seq.PrintFold()
        else:
            results.loc[id,'opt_score'] = score 
            results.loc[id,'c_score'] = "The folding %s has wrong length." % (fold)
    
    results.to_excel('results.xlsx')

     
if __name__ == "__main__":
    main()
    
    