import random
import sys
import math
from numpy.random import choice

nucleotides = ["A", "T", "C", "G"]

def create_binding_sites(sequence_count, motif_matrix):
    binding_sites = []
    for i in range(sequence_count):
        seq = []
        for j in range(len(motif_matrix[0])):
             seq.append(choice(nucleotides, p=[motif_matrix[k][j] for k in range(len(motif_matrix))]))
        binding_sites.append("".join(seq))
    return binding_sites

def create_random_motif(icpc, length):
    """ Creates a 4 x `length` matrix, where the sum of each column of the matrix is `icpc` """
    motif = [[0 for cols in range(length)] for rows in range(4)]
    for i in range(4):
        for j in range(length):
            if i == 3:
                motif[i][j] = icpc - sum(row[j] for row in motif)
            else:
                motif[i][j] = random.uniform(0, (icpc - sum(row[j] for row in motif)))

    return motif

def create_sequences(count, length):
    """ Creates `count` sequences of `length` each. Returns an array of length `count`"""
    seqs = []
    for i in range(count):
        seq = [random.choice(nucleotides) for j in range(length)]
	seqs.append("".join(seq))
    return seqs

def main():
    if len(sys.argv) <= 4:
        print "You must provide four numbers for the benchmark"
        return

    icpc = float(sys.argv[1])
    ml = int(sys.argv[2])
    sl = int(sys.argv[3])
    sc = int(sys.argv[4])

    sequences = create_sequences(sc, sl)
    motif = create_random_motif(icpc, ml)
<<<<<<< HEAD
    binding_sites = create_binding_sites(sc, motif)
    print(binding_sites)

=======

    print (motif)

>>>>>>> 4ae6a05... Fix summation error
if __name__ == "__main__":
    main()
