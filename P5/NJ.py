import numpy as np
from parsing import *
import sys

def calculateN(d,nodes):
    s = len(d)
    N = np.zeros((s,s))
    
    for i in range(s):
        for j in range(s):
            if i != j:
                sumI = (1/ (s-2)) * sum(d[i,label] for label in nodes.values())
                sumJ = (1/ (s-2)) * sum(d[j,label] for label in nodes.values())
                N[i][j] = d[i,j] - (sumI + sumJ)
            else:
                N[i][j] = 0 
    return N


def minimum_entry(N):
    return np.unravel_index(N.argmin(), N.shape)


def getNewEdges(i,j,d):
    
    sumI = sum(d[i,k] for k in range(len(d)))
    sumJ = sum(d[j,k] for k in range(len(d)))
    
    dki = round(0.5 * (d[i,j] + sumI  - sumJ),3)
    dkj = round(0.5 * (d[i,j] + sumJ  - sumI),3) #- dki 
    return (dki,dkj)


def updateDistanceMatrix(a,b,m,keys):
    ## get the indices of the nodes to merge
    i = keys[a]
    j = keys[b]
    
    ## update the new keys
    updated_nodes = keys.copy()
    ### remove a and b
    updated_nodes.pop(a)
    updated_nodes.pop(b)
    ### update the indices
    updated_nodes = {k: i for i, k in enumerate(sorted(updated_nodes.keys(), key=lambda x: updated_nodes[x]))}
    
    
    size_nd = len(m)-1
    nd = np.zeros((size_nd,size_nd))
    ## update distance matrix for the nodes that are not involved in the merging operation, keep the distances
    for node, upd_index in updated_nodes.items():
        index_in_m = keys.get(node)
        for node_ , upd_index_ in updated_nodes.items():
            index_in_m_ = keys.get(node_)
            nd[upd_index,upd_index_] = m[index_in_m,index_in_m_]
            
    ### generate a new index for the new node
    updated_nodes[a+b] = len(updated_nodes)    
    
    ## compute the distance from all the nodes to the merged node
    for node, upd_index in updated_nodes.items():
        # get the index in the original distance matrix 
        index_in_m = keys.get(node)
        # if the node exists in the matrix (it is not the merged node), get the distance
        if index_in_m is not None:
            nd[-1][upd_index] = (m[i][index_in_m] + m[j][index_in_m] - m[i][j]) / 2.
            nd[upd_index][-1] =(m[i][index_in_m] + m[j][index_in_m] - m[i][j]) / 2.

    return nd,updated_nodes,a+b
    


def NeighbourJoining(d,nodes):
    
    S = len(d)
    T = {}
    # d -> np array + dictionary of nodes and indices
    while S > 3:
        ## 1.a Compute N
        N = calculateN(d,nodes)
        ## 1.b Find the min in N
        lowestPair = minimum_entry(N)

        ## 2. Add a new node k to the tree T
        ## 3. add edges with weights
        i = lowestPair[0]
        j = lowestPair[1]
        newEdges = getNewEdges(i,j,d)
        #nodes[(i,j)] = newEdges
        
        node_a = list(nodes.keys())[i]
        node_b = list(nodes.keys())[j]
        ## 4. Update the ds matrix by deleting rows and columns corresponding
        ## to i and j and adding a new row and column for the new taxon k
        d,nodes,new_node = updateDistanceMatrix(node_a,node_b,d,nodes)
        
        ## save Newick format
        ## print("Merging: ({}:{},{}:{})".format(node_a, newEdges[0], node_b, newEdges[1]))
        T[new_node] = f"({node_a}:{newEdges[0]},{node_b}:{newEdges[1]})"
        
        S-=1
    
    ## termination
    i, j, m = 0,1,2
    remaining_nodes = list(nodes.keys())
    v_i = round(0.5 * (d[i,j] + d[i,m] - d[j,m]),3)
    v_j = round(0.5 * (d[i,j] + d[j,m] - d[i,m]),3)
    v_m = round(0.5 * (d[i,m] + d[j,m] - d[i,j]),3)
    
    newick = f"({remaining_nodes[i]}:{v_i},{remaining_nodes[j]}:{v_j},{remaining_nodes[m]}:{v_m})"+';'
    
    for clade, newick_format in T.items():
        newick = newick.replace(clade,newick_format)
    
    return newick
    


def NJ(phy_file, outputfile = None):
    
    D, nodes = parse_phy_file(phy_file)
  
    tree = NeighbourJoining(D,nodes)
    

    with open(outputfile, "w") as f:
        f.write(tree)
        
    ## return tree

def main():
    
    if len(sys.argv) == 3:
        NJ(sys.argv[1],sys.argv[2])
    else: 
        outputfile = "output/"+ sys.argv[1][-10:-4] + ".nwk"
        NJ(sys.argv[1],outputfile)
        

if __name__ == "__main__":
    main()
   
    # phy_file = "P5/unique_distance_matrices/89_Adeno_E3_CR1.phy"
    # outputfile  = "P5/output/89_Adeno_E3_CR1.nwk"
    # print(NJ(phy_file,outputfile))
