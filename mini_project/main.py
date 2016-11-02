import random
import sys

nucleotides = ["A", "T", "C", "G"]

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

if __name__ == "__main__":
    main()
