import random
import sys
import math
import numpy as np
import os
import copy
import timeit
import multiprocessing

from numpy.random import choice

nucleotides = {"A" : 0, "T" : 1, "C" : 2, "G" : 3}

set_count = 1

def calculate_information_content(pwm):
    ic = 0
    for j in xrange(len(pwm[0])):
        for i in range(len(pwm)):
            ic += pwm[i][j] * math.log(((pwm[i][j]+.00001)/.25), 2)
    return ic 

def output_elapsed_time(elapsed):
    """
    Output the elapsed time
    """
    with open("elapsed_time.txt", "w+") as output_file:
        output_file.write(str(elapsed))

def output_predicted_motif(motif, ml):
    """
    Writes the motif matrix out to predicted_motif.txt given the specific format
    """
    rows = len(motif)
    cols = len(motif[0])
    with open("predictedmotif.txt", "w+") as output_file:
        output_file.write(">PREDICTED-MOTIF" + "\t" + str(ml))
        for i in xrange(rows):
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

def create_pwm(sequences, icpc):
    pwm = [[], [], [], []]
    for i in xrange(len(sequences[0])):
        nucleotide_count = [0, 0, 0, 0]
        for j in sequences:
            nucleotide_count[nucleotides[j[i]]] += 1.0
        for t in xrange(4):
            pwm[t].append(nucleotide_count[t])
    for i in xrange(len(pwm[0])):
        total = sum([row[i] for row in pwm])
        for t in xrange(4):
            pwm[t][i] /= total
            pwm[t][i] *= icpc
    return pwm

def check_score(sequence, pwm):
    score = 0
    for i in xrange(len(sequence)):
        score += (pwm[nucleotides[sequence[i]]][i] + 0.00001)/.25     # This calculates Qx/Px, which the score
    return score

def find_best_location(pwm, sequence, ml):
    score_list = []
    for i in xrange(len(sequence) - ml + 1):
        score = check_score(sequence[i:i+ml], pwm)
        score_list.append(score)
    total = sum(score_list)
    score_list = map(lambda x: x/total, score_list)     # Convert the score list into a probability distribution by dividing each element by the total
    x = np.random.choice(np.arange(0, len(score_list)), p=score_list)   # Sample from the probability distribution
    # max_score = max(score_list)
    # max_list = [i for i, j in enumerate(score_list) if j == max_score]
    return x

def create_motif(sequences, ml, icpc):
    """Finds motif of length `ml` in `sequences` using Gibbs Motif Sampling"""
    max_ic = 0
    current_best_sites = []
    for x in range(3):     #Go through the process multiple times
        prev_sites = []
        sites = [random.randint(0, len(i) - ml) for i in sequences]
        while(prev_sites != sites):
            prev_sites = copy.deepcopy(sites)
            i = random.randint(0, len(prev_sites) - 1)
            # for i in range(len(prev_sites)):
            part = [sequences[t][sites[t]:sites[t] + ml] for t in xrange(len(sites)) if t != i]    #Get the relevant sections of the sequences for the pwm for all sequence other than i
            pwm = create_pwm(part, icpc)
            best_pos = find_best_location(pwm, sequences[i], ml)
            sites[i] = best_pos
            
            full_motifs = [sequences[t][sites[t]:sites[t] + ml] for t in xrange(len(sites))]
            full_pwm = create_pwm(full_motifs, icpc)
            ic = calculate_information_content(full_pwm)

            if(ic > max_ic):    #Pick the sites with the Highest ic content according to lecture
                max_ic = ic
                current_best_sites = copy.deepcopy(sites)
    return create_pwm([sequences[t][current_best_sites[t]:current_best_sites[t] + ml] for t in xrange(len(current_best_sites))], icpc), current_best_sites

def find_motifs(items):
    icpc = items[0]
    set_count = items[1]

    print("Set: %d" % set_count)
    base_path = os.getcwd()
    os.chdir("data/set%d" % set_count)

    ml = read_motif_length()
    sequences = read_sequences()

    start_time = timeit.default_timer()
    motif, plant_sites = create_motif(sequences, ml, icpc)

    # Transpose according to grading rubric
    motif = np.transpose(unnormalize_motif(motif, len(sequences)))
    elapsed = timeit.default_timer() - start_time
    output_predicted_plant_sites(plant_sites)
    output_predicted_motif(motif, ml)
    output_elapsed_time(elapsed)

    os.chdir(base_path)

def main():
    global set_count
    if not os.path.exists("data"):
        print "Exiting the program. Data benchmarks have not been generated" 
        return
    else:
        default_icpc = 2
        default_ml = 8
        default_sc = 10

        icpc_vals = [1, 1.5]
        motif_lengths = [6, 7]
        sequence_counts = [5, 20]

        args_queue = []

        for i in range(1, 11):
            args_queue.append([default_icpc, set_count])
            set_count += 1

        for icpc in icpc_vals:
            for i in range(1, 11):
                args_queue.append([icpc, set_count])
                set_count += 1

        for ml in motif_lengths:
            for i in range(1, 11):
                args_queue.append([default_icpc, set_count])
                set_count += 1

        for sc in sequence_counts:
            for i in range(1, 11):
                args_queue.append([default_icpc, set_count])
                set_count += 1

        worker_pool = multiprocessing.Pool(processes=2)
        worker_pool.map(find_motifs, args_queue)
        
if __name__ == "__main__":
    main()
