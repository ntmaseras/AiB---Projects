import numpy as np
import itertools
from parsing import *

def calculateN(d):
    s = len(d)
    N = np.zeros((s,s))
    #d,_ = getValues(d)
    labels = list(d.keys())
    for i in range(s):
        for j in range(s):
            if i != j:
                sumI = (1/ (s-2)) * sum(d[labels[i]][label] for label in labels)
                sumJ = (1/ (s-2)) * sum(d[labels[j]][label] for label in labels)
                N[i][j] = d[labels[i]][labels[j]] - (sumI + sumJ)
            else:
                N[i][j] = 0 
    return N


def minimum_entry(N):
    return np.unravel_index(N.argmin(), N.shape)


def getNewEdges(i,j,d):
    d,keys = getValues(d)
    sumI = 0
    sumJ = 0
    for k in range(len(d)):
        sumI += d[keys[i]][k]
        sumJ += d[keys[j]][k]

    dki = 0.5 * (d[keys[i]][keys[j]] + sumI  - sumJ)
    dkj = 0.5 *(d[keys[i]][keys[j]] + sumJ - sumI) #- dki 
    return (dki,dkj)

def updateDistanceMatrix(a,b,D,nodes):
     
    D,keys = getValues(D)
    original_set_of_nodes = keys.copy()
    original_set_of_nodes.pop(a)
    original_set_of_nodes.pop(b)
    
    m = np.vstack(D)
    i = keys[a]
    j = keys[b]
    
    
    size_nd = len(m)-1
    nd = np.empty((size_nd,size_nd))
    nd[:] = np.NaN
    

    for v in range(0,len(nd)-1):
        if v != i or v != j:
            for w in range(0,len(nd)-1):
                if w != i or  w != j:
                    nd[v,w] = m[v,w]
      
    for node in original_set_of_nodes.values():
        for d in range(0,len(nd)):
            if len(nd)-1 == d:
                nd[-1][-1] = 0
            else:
                nd[-1][d] = (m[i][node] + m[j][node] - m[i][j]) / 2.
                nd[d][-1] =(m[i][node] + m[j][node] - m[i][j]) / 2.

                
   
    
    return nd
    
def NeighbourJoining(d):
    labels = list(d.keys())
    nodes = {}
    S = len(d)
    T = []
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
        T.append("({}:{},{}:{})".format(i, newEdges[0], j, newEdges[1])+",")
        ## 4. Update the ds matrix by deleting rows and columns corresponding
        ## to i and j and adding a new row and column for the new taxon k
        d = updateDistanceMatrix(i,j,d,list(nodes.values()))
       
        S-=1
        
        return T
    

if __name__ == "__main__":
    D, _ = parse_matrix_and_gap("P5/example_slide4.phy")
    NeighbourJoining(D)
