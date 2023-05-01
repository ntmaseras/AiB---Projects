import random
import numpy as np

def initialize_lattice(size):
    """
    Initializes a 2D lattice of a given size with all elements set to 0.
    """
    lattice = np.zeros((size, size), dtype=int)
    return lattice

def initialize_protein(hp_string):
    """
    Initializes a protein chain from an HP string, placing the first amino acid
    at the center of the lattice and the remaining amino acids in a diagonal line.
    """
    size = len(hp_string) * 2 + 1
    lattice = initialize_lattice(size)
    x = size // 2
    y = size // 2
    for amino_acid in hp_string:
        if amino_acid == "H":
            lattice[x][y] = 1
        x += 1
        y += 1
    return lattice

def score_protein(protein):
    """
    Computes the energy of a folded protein using the HP model scoring function.
    """
    energy = 0
    size = protein.shape[0]
    for x in range(size):
        for y in range(size):
            if protein[x][y] == 1:
                if x > 0 and protein[x-1][y] == 1:
                    energy -= 1
                if y > 0 and protein[x][y-1] == 1:
                    energy -= 1
                if x < size-1 and protein[x+1][y] == 1:
                    energy -= 1
                if y < size-1 and protein[x][y+1] == 1:
                    energy -= 1
    return energy

def fold_protein(protein, iterations):
    """
    Implements the 1/4 approximation algorithm for folding a protein chain in a given
    number of iterations.
    """
    size = protein.shape[0]
    for i in range(iterations):
        x = random.randint(1, size - 2)
        y = random.randint(1, size - 2)
        amino_acid = protein[x][y]
        neighbor_sum = (protein[x-1][y] + protein[x+1][y] +
                        protein[x][y-1] + protein[x][y+1])
        if amino_acid == 0 and neighbor_sum == 1:
            protein[x][y] = 1
        elif amino_acid == 1 and neighbor_sum in (0, 2):
            protein[x][y] = 0
    return protein
def visualize_protein(protein):
    """
    Generates an ASCII art representation of a folded protein chain.
    """
    protein_str = ""
    for row in protein:
        protein_str += "".join(["H" if val == 1 else "." for val in row])
        protein_str += "\n"
    return protein_str

# Example usage:
hp_string = "HPHPPHHPHPPH"
protein = initialize_protein(hp_string)
folded_protein = fold_protein(protein, 10000)
energy = score_protein(folded_protein)
print("Folded protein:")
print(visualize_protein(folded_protein))
print(f"Energy: {energy}")