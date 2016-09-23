import sys
import time
import random
import subprocess

nucleotides = ["A", "C", "G", "T"]

def create_sequences(length):
    s1 = [random.choice(nucleotides) for i in range(length)]
    s2 = list(s1)
    s3 = list(s1)
    for i in range(length/10):
        cur = s2
        for j in range(2):
            pos = random.randint(0, len(cur) - 1)
            change = random.randint(0, 1)
            if change:
                filtered = [n for n in nucleotides if n != cur[pos]]
                cur[pos] = random.choice(filtered)
            else:
                del cur[pos]
            cur = s3
            
    return "".join(s2), "".join(s3)
    
def save_sequence_as_fasta(seq, file_name):
    with open(file_name, 'w') as fasta_file:
        fasta_file.write(">Sequence %s\n" % (file_name))
        fasta_file.write(seq)

def run_align(seq1_file, seq2_file, subs_file, gap_penalty):
    p = subprocess.Popen(["python", "align.py", seq1_file, seq2_file, subs_file, gap_penalty], stdout=subprocess.PIPE)
    output = p.communicate()
    arr = output[0].split("\n")
    for i in arr:
        if i != "":
            print(i)

def main():
    if len(sys.argv) != 2:
        print("You must provide a length for your sequence")
    else:
        random.seed(time.time())
        length = sys.argv[1]
        seq1, seq2 = create_sequences(int(length))
        save_sequence_as_fasta(seq1, 'seq1.fasta')
        save_sequence_as_fasta(seq2, 'seq2.fasta')
        run_align('seq1.fasta', 'seq2.fasta', 'subs.txt', '-500')


if __name__ == "__main__":
    main()
