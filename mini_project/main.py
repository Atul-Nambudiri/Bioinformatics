import random
import sys
import math
from numpy.random import choice

nucleotides = ["A", "T", "C", "G"]

def output_plant_sites(plant_sites):
    """ Outputs the `plant_sites` array to sites.txt, one site per line """ 
    with open("sites.txt", "w+") as output_file:
        output_file.write("\n".join(str(i) for i in plant_sites))

def output_sequences(sequences):
    """ Prints the `sequences` out to "sequences.fa", in FASTA Format """
    with open("sequences.fa", "w+") as output_file:
        output_file.write(">" + "\n>".join(sequences))

def create_binding_sites(sequence_count, motif_matrix, length):
    """ Create `sequence_count` `binding-sites` of length `length`, by sampling the `motif-matrix` """
    binding_sites = []
    for i in range(sequence_count):
        seq = []
        for j in range(length):
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

def plant_site(seqs, sites):
    """ Over-rides a substring from each sequence in `seqs` with a seq from `sites` """
    new_seqs = []
    plant_sites = []
    for i in range(len(seqs)):
        idx = random.randint(0, len(seqs[i]) - len(sites[i]) - 1)
        new_seqs.append(seqs[i][:idx] + sites[i] + seqs[i][idx+len(sites[i]):])
        plant_sites.append(idx)
    return new_seqs, plant_sites

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
    binding_sites = create_binding_sites(sc, motif, ml)
    new_seqs, plant_sites = plant_site(sequences, binding_sites)

    print(sequences)
    print(binding_sites)
    print(new_seqs)
    print(plant_sites)

    output_sequences(new_seqs)
    output_plant_sites(plant_sites)




if __name__ == "__main__":
    main()
