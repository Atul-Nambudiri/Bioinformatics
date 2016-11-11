import random
import sys
import math
import numpy as np
import os

from numpy.random import choice

nucleotides = ["A", "T", "C", "G"]

def output_predicted_motif(motif, ml, idx):
    """
    Writes the motif matrix out to predicted_motif.txt given the specific format
    The `idx` parameter is to be used when creating the benchmark datasets
    """
    rows = len(motif)
    cols = len(motif[0])
    with open("predictedmotif.txt", "w+") as output_file:
        output_file.write(">MOTIF" + str(idx) + "\t" + str(ml))
        for i in range(rows):
            output_file.write("\n")
            for j in range(cols):
                output_file.write(str(motif[i][j]) + "\t")
        output_file.write("\n<")

def output_predicted_plant_sites(plant_sites):
    """ Outputs the `plant_sites` array to predictedsites.txt, one site per line """
    with open("predictedsites.txt", "w+") as output_file:
        output_file.write("\n".join(str(i) for i in plant_sites))

def read_motif_length():
    """ Reads and returns the motif length from motiflength.txt """
    ret = -1
    with open("motiflength.txt", "r") as output_file:
        ret = int(output_file.read())
    return ret

def read_sequences():
    """ Reads the `sequences` out to "sequences.fa """
    sequences = []
    with open("sequences.fa", "r") as output_file:
        for i in output_file:
            sequences.append(i[1:len(i)-1])
    return sequences


def find_motifs():
    motif_length = read_motif_length()
    sequences = read_sequences()
    print(motif_length)
    print(sequences)

def main():
    if not os.path.exists("data"):
        print "Exiting the program. Data benchmarks have not been generated" 
        return
    else:
        if len(sys.argv) < 2:
            print "You must provide the name of the benchmark to test"
            return
        benchmark = sys.argv[1]
        base_path = os.getcwd()
        os.chdir("data/%s" % benchmark)
        find_motifs()
        os.chdir(base_path)

if __name__ == "__main__":
    main()
