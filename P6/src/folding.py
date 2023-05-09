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
    
    evens_b = True
    ## get the max 
    if len(matching_evens) > len(matching_odds):
        matching = matching_evens
        evens = list(matching_evens.keys())
        odds = list(matching_evens.values())
        i,j = i_even, j_even
    else:
        matching = matching_odds
        evens = list(matching_odds.values())
        odds = list(matching_odds.keys())
        evens_b = False
        
    p = (i+j) // 2
    
    
    return dict(sorted(matching.items())),evens, odds, p




def fold_protein(evens, odds, p, n):
    
    """
    Implements the 1/4 approximation algorithm for folding a protein chain in a given
    number of iterations.
    """
    
    fold = (n-1) * ['F']
    fold[p] = 'r'
    fold[p+1] = 'r'
    ## ALL PROBLEMS ARE AROUND P
    for i in range(0, len(evens) - 1):
        current =  evens[i]
        next = evens[i + 1]
        diff = next - current
        if diff > 3:   
            fold[current] = 'L'
            roof = int(current + diff // 2)
            fold[roof - 1] =  'R'
            fold[roof] =  'R'
            fold[next - 1] = 'L'
            if current < p < next:
                ## if # of bases before the split is > 3, build the house
                if (p - current) >= 3: 
                    fold[next - 1] = 'F'
                elif (p - current) < 3:
                    fold[current] = 'f'
                    fold[p] = 'r'
                    fold[p+1] = 'f'
                elif current == p + 1:
                    fold[current] = 'F'

    for i in range(0, len(odds) - 1):
        current =  odds[i]
        next = odds[i + 1]
        diff = next - current

        if diff > 3:       
            fold[current] = 'L'
            roof = int(current + diff // 2)
            fold[roof - 1] =  'R'
            fold[roof] =  'R'
            fold[next - 1] = 'L'
            if current < p < next:
                ## if # of bases before the split is > 3, build the house
                if (p - current) >= 3: 
                    fold[next - 1] = 'F'
                elif (p - current) < 3:
                    fold[current] = 'f'
                    fold[p] = 'r'
                    fold[p+1] = 'f'
                elif current == p + 1:
                    fold[current] = 'F'

                #elif p == next - 1:
            

    return ''.join(fold)



def main():

    filename = 'P6/src/sequences.txt'
    sequences = {}
    with open(filename, "r") as f:
        for i, line in enumerate(f.readlines()):
            seq = line.split(' ')[0]
            score = int(line.split(' ')[1].replace('\n',''))
            sequences[i] = [seq,score]
    
    results = pd.DataFrame()
    for id, [seq,score] in sequences.items():
       
        matching,evens, odds, p = initialize_protein(seq.lower())
        print(matching)
        print(len(seq))
        fold = fold_protein(evens, odds, p, len(seq)) 
        print(id ,fold, p)
        
        results.loc[id,'seq'] = seq 
        seq = HPFold(seq)
        absfold = make_absfold(fold)
        relfold = make_relfold(fold)
        if len(absfold) == len(seq) - 1:
            seq.SetAbsFold(absfold)
        elif len(relfold) == len(seq) - 1:
            seq.SetRelFold(relfold)
        
        results.loc[id,'opt_score'] = score 
        results.loc[id,'aprox_score'] = -len(matching)
        score, illegal = seq.PrintFold()
        results.loc[id,'capi'] = -score
        #results.loc[id,'1/4 opt'] = 0.25*score 

    
        results.loc[id,'Illegal?'] = illegal

    
    results.to_excel('results.xlsx')

     
if __name__ == "__main__":
    main()
    
    