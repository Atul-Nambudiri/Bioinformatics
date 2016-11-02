import sys

def read_fasta_file(file_name):
    sequence_string = ""
    with open(file_name, 'r') as fasta_file:
        for line in fasta_file:
            if line[0] != ">":
                res = line.replace("\n", "").replace("\r", "")
                sequence_string += res
    return sequence_string

def read_substitution_matrix_file(file_name):
    matrix = {}
    with open(file_name, 'r') as subs:
        for i, line in enumerate(subs):
            res_arr = line.split()
            if i == 0:
                for j in range(0, 4):
                    matrix[res_arr[j]] = {}
            else:
                matrix[res_arr[0]]['A'] = int(res_arr[1])
                matrix[res_arr[0]]['C'] = int(res_arr[2])
                matrix[res_arr[0]]['G'] = int(res_arr[3])
                matrix[res_arr[0]]['T'] = int(res_arr[4])
    return matrix

def init_matrix(seq1, seq2, gap):
    matrix = []
    for i in range(len(seq1) + 2):
    	matrix.append([0 for j in range(len(seq2) + 2)])
    matrix[1][1] = 0  # Initialize matrix
    for i in range(len(matrix)):
	if i == 0:
            continue
	if i == 1:
            matrix[i][0] = "s1"
	else:
            matrix[i][0] = seq1[i -2]
            matrix[i][1] = matrix[i -1][1] + gap
    for j in range(len(matrix[0])):
	if j == 0:
            continue
	if j == 1:
            matrix[0][j] = "s2"
        else:
            matrix[0][j] = seq2[j -2]
            matrix[1][j] = matrix[1][j-1] + gap
    return matrix

def fill_matrix(seq1, seq2, sub_matrix, gap):
    matrix = init_matrix(seq1, seq2, gap)
    for i in range(2, len(seq1)+2):
        for j in range(2, len(seq2)+2):
            diag = matrix[i-1][j-1] + (sub_matrix[seq1[i-2]][seq2[j-2]]) 
            top = matrix[i-1][j] + gap
            left = matrix[i][j-1] + gap
            matrix[i][j] = max([diag, top, left])
    return matrix

def print_aligned_sequences(seq1, seq2, solved_matrix, sub_matrix, gap):
    as1 = ""
    as2 = ""
    i = len(seq1)
    j = len(seq2)
    while(i > 0 and j > 0):
        if seq1[i-1] == seq2[j-1] and solved_matrix[i][j] == solved_matrix[i+1][j+1] - (sub_matrix[seq1[i-1]][seq2[j-1]]):
            as1 = seq1[i - 1] + as1
            as2 = seq2[j -1] + as2
            i -= 1
            j -= 1
        elif seq1[i-1] != seq2[j-1] and solved_matrix[i][j] == solved_matrix[i+1][j+1] - (sub_matrix[seq1[i-1]][seq2[j-1]]):
            as1 = seq1[i - 1] + as1
            as2 = seq2[j -1] + as2
            i -= 1
            j -= 1
        elif solved_matrix[i+1][j] == solved_matrix[i+1][j+1] - gap:
            as1 = "-" + as1
            as2 = seq2[j -1] + as2
            j -= 1 
        elif solved_matrix[i][j+1] == solved_matrix[i+1][j+1] - gap:
            as1 = seq1[i - 1] + as1
            as2 = "-" + as2
            i -= 1

    while(i > 0):
        as1 = seq1[i -1] + as1
        as2 = "-" + as2
        i -= 1
    while(j > 0) :
        as1 = "-" + as1
        as2 = seq2[j-1] + as2
        j -= 1
   
    score = 0 

    for i in range(len(as1)):
        if as1[i] == "-" or as2[i] == "-":
            score += gap
        else:
            score += sub_matrix[as1[i]][as2[i]]
    
    print("The optimal alignment between given sequences has score %d" % (score))
    print(as1)
    print(as2)


def print_matrix(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print '\n'.join(table)

def main():
    if len(sys.argv) != 5:
        print("You must pass in two fasta files, a substitution matrix file, and a negative gap penalty to this program")
    else:
        sub_matrix = read_substitution_matrix_file(sys.argv[3]) 
        seq1 = read_fasta_file(sys.argv[1])
        seq2 = read_fasta_file(sys.argv[2])
        matrix = fill_matrix(seq1, seq2, sub_matrix, int(sys.argv[4]))
        print_aligned_sequences(seq1, seq2, matrix, sub_matrix, int(sys.argv[4]))

if __name__ == "__main__":
    main()
		
		
		
