from Bio import Phylo
from Bio.Phylo.Consensus import _BitString
import sys 
# Parse the Newick tree
def parse_newick(tree): #reads file and returns tree
    return Phylo.read(tree, "newick")

def get_subtrees(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString("".join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
        
    return bitstrs

def get_rf_distance(tree1, tree2):
    subtrees1 = get_subtrees(tree1)
    subtrees2= get_subtrees(tree2)
    return len(subtrees1.union(subtrees2)) - len(subtrees1.intersection(subtrees2))


def main():
    ## test
    tree1 = parse_newick('Testdata/tree1.new')
    tree2 = parse_newick('Testdata/tree2.new')
    ## dist = get_rf_distance(tree1,tree2)
    ## print("Testdata distance: ",dist)
    ####
    #tree1 = parse_newick(sys.argv[1])
    #tree2 = parse_newick(sys.argv[2])
    dist = get_rf_distance(tree1,tree2)
    print("Given trees distance: ",dist)
    
    
if __name__ == '__main__':
    main()