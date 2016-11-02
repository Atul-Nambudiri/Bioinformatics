import random
import sys

nucleotides = ["A", "T", "C", "G"]

def create_sequences(count, length):
    """ Creates `count` sequences of `length` each. Returns an array of length `count`"""
    seqs = []
    for i in range(count):
        seq = ""
        for j in range(length):
            seq += random.choice(nucleotides)
        seqs.append(seq)

    return seqs

def main():
    if len(sys.argv) < 4:
        print "You must provide four numbers for the benchmark"
        return

    icpc = int(sys.argv[1])
    ml = int(sys.argv[2])
    sl = int(sys.argv[3])
    sc = int(sys.argv[4])

    sequences = create_sequences(sc, sl)
    print (sequences)


if __name__ == "__main__":
    main()
