import random
import sys
import math
import numpy as np

from numpy.random import choice

nucleotides = ["A", "T", "C", "G"]

def output_motif(motif, ml, idx):
    """
    Writes the motif matrix out to motif.txt given the specific format
    The `idx` parameter is to be used when creating the benchmark datasets
    """
    rows = len(motif)
    cols = len(motif[0])
    with open("motif.txt", "w+") as output_file:
        output_file.write(">MOTIF" + str(idx) + "\t" + str(ml))
        for i in range(rows):
            output_file.write("\n")
            for j in range(cols):
                output_file.write(str(motif[i][j]) + "\t")
        output_file.write("\n<")

def unnormalize_motif(motif, sc):
    """ Returns the unnormalized matrix motif """
    rows = len(motif)
    cols = len(motif[0])
    new_mat = [[0 for i in range(cols)] for j in range(rows)]
    for i in range(rows):
        for j in range(cols):
            new_mat[i][j] = int(round(motif[i][j] * sc))

    return new_mat

def output_motif_length(ml):
    """ Outputs the motif length `ml` to a file motiflength.txt """
    with open("motiflength.txt", "w+") as output_file:
        output_file.write(str(ml) + "\n")

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

    # Transpose according to grading rubric
    unnormalized_motif = np.transpose(unnormalize_motif(motif, sc))

    output_sequences(new_seqs)
    output_plant_sites(plant_sites)
    output_motif(unnormalized_motif, ml, 1)
    output_motif_length(ml)

if __name__ == "__main__":
    main()
