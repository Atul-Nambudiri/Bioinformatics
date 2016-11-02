import math 
from numpy.random import choice

nucleotides = ["A", "T", "C", "G"]

def main():
    sequence_count = 20
    matrix = [[.8, .05, .1],
              [.1, .1, .7],
              [.05, .05, .1],
              [.05, .8, .1]]
    binding_sites = []
    for i in range(sequence_count):
        seq = []
        for j in range(len(matrix[0])):
             seq.append(choice(nucleotides, p=[matrix[k][j] for k in range(len(matrix))]))
        binding_sites.append("".join(seq))




if __name__ == "__main__":
    main()