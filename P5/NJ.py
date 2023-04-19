import sys
import numpy as np
import scipy as scipy
import itertools
from parsing import *

def calculateN(d):
    s = len(d)
    N = np.zeros((s,s))
    labels = list(d.keys())
    for i in range(s):
        for j in range(s):
            sumI = (1/ (s-2)) * sum(d[labels[i]][label] for label in labels if label != labels[i])
            sumJ = (1/ (s-2)) * sum(d[labels[j]][label] for label in labels if label != labels[j])
            N[i][j] = d[labels[i]][labels[j]] - (sumI + sumJ)
    return N


def minimum_entry(N):
    min_val = np.inf
    min_index = None
    for i, j in itertools.combinations(range(len(N)), 2):
        if N[i][j] < min_val:
            min_val = N[i][j]
            min_index = (i, j)
    return min_index


def getNewEdges(i,j,d):
    d,keys = getValues(d)
    sumI = 0
    sumJ = 0
    for k in range(len(d)):
        sumI += d[keys[i]][k]
        sumJ += d[keys[j]][k]

    dki = 0.5 * (d[keys[i]][keys[j]]) - (sumI + sumJ)
    dkj = 0.5 * (d[keys[i]][keys[j]]) - dki 
    return (dki,dkj)

def updateDistanceMatrix(i,j,d):
    d,keys = getValues(d)
    r = len(d)
    nd = 0
    return nd
    
def NeighbourJoining(d):
    labels = list(d.keys())
    nodes = {}
    S = len(d)
    while S > 3:
        ## 1.a Compute N
        N = calculateN(d)
        ## 1.b Find the min in N
        lowestPair = minimum_entry(N)

        ## 2. Add a new node k to the tree T
        ## 3. add edges with weights
        i = labels[lowestPair[0]]
        j = labels[lowestPair[1]]
        newEdges = getNewEdges(i,j,d)
        nodes[(i,j)] = newEdges
        print("Merging: ({}:{},{}:{})".format(i, newEdges[0], j, newEdges[1]))

        ## 4. Update the d matrix by deleting rows and columns corresponding
        ## to i and j and adding a new row and column for the new taxon k
        a = updateDistanceMatrix(i,j,d)
       
        S-=1
        
    

if __name__ == "__main__":
    D, gap_cost = parse_matrix_and_gap("P5/test.phy")
    NeighbourJoining(D)
