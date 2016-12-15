import sys
import os
from math import log
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


set_count = 1

"""
This file contains all of the evaluation logic for the final project
"""

def compute_overlapping_sites(predicted, actual):
    overlap = 0
    for site in predicted:
        if site in actual:
            overlap += 1
    return overlap

def compute_overlapping_pos(predicted, actual, ml):
    p = []
    a = []

    for site in predicted:
        for pos in range(site+1, site+ml+1):
            p.append(pos)

    for site in actual:
        for pos in range(site+1, site+ml+1):
            a.append(pos)

    overlap = 0
    for pos in p:
        if pos in a:
            overlap +=1
    return overlap

def parse_sites(file_name):
    sites = []
    with open(file_name, "r") as f:
        for line in f:
            sites.append(int(line))
    return sites

def parse_ml(file_name):
    ml = 0
    with open(file_name, "r") as f:
        ml = int(f.readline())
    return ml


def get_path(folder, f):
    return "/".join([folder, f])

def parse_motif(motif_file):
    """ Reads and returns the motif from `motif_file` """
    motif = []
    with open(motif_file, "r") as input_file:
        for i in input_file:
            if i[0] not in ["<", ">"]:
                motif.append([float(t) + .00001 for t in (i.replace("\t\n", "").split("\t"))]) #Add a small amount to each value to prevent division by 0 errors
    return motif

def calculate_position_relative_entropy(actual_position, predicted_position):
    """ Calculates the relative entropy between the two motifs
        at one position `actual_position`, `predicted_position`
        Based on page 6 of this: http://www.seas.upenn.edu/~cis535/Lab/motiffinding1.pdf
    """
    actual_positon_total_score = float(sum(actual_position))
    predicted_position_total_score = float(sum(predicted_position))
    entropy = 0
    for i in xrange(len(actual_position)):
        actual_normalized = actual_position[i]/actual_positon_total_score
        predicted_normalized = predicted_position[i]/predicted_position_total_score
        entropy += (predicted_normalized * log(predicted_normalized/actual_normalized, 2))
    return entropy

def calculate_relative_entropy(actual_motif, predicted_motif):
    entropy = 0
    for i in xrange(len(actual_motif)):
        entropy += calculate_position_relative_entropy(actual_motif[i], predicted_motif[i])
    return entropy

def benchmark_helper(icpc, ml, sc):
    global set_count

    motif_file = "motif.txt"
    predicted_motif_file = "predictedmotif.txt"

    ml_file = "motiflength.txt"

    sites_file = "sites.txt"
    predicted_sites_file = "predictedsites.txt"

    base_path = os.getcwd()
    os.chdir("data/set%d" % set_count)

    actual_motif = parse_motif(motif_file)
    predicted_motif = parse_motif(predicted_motif_file)

    ml = parse_ml(ml_file)
    sites = parse_sites(sites_file)
    p_sites = parse_sites(predicted_sites_file)

    entropy = calculate_relative_entropy(actual_motif, predicted_motif)
    overlap_sites = compute_overlapping_sites(p_sites, sites)
    overlap_pos = compute_overlapping_pos(p_sites, sites, ml)

    with open("entropy.txt", "w+") as output_file:
        output_file.write(str(entropy))

    with open("overlap_pos.txt", "w+") as output_file:
        output_file.write(str(overlap_pos))

    with open("overlap_sites.txt", "w+") as output_file:
        output_file.write(str(overlap_sites))

    #print "======== DATA FOR SET%d =========" % set_count
    #print("Entropy: %f:" % entropy)
    #print("Number of positional overlaps: %d" % overlap_pos)
    #print("Number of site overlaps: %d" % overlap_sites)
    #print ""

    os.chdir(base_path)
    set_count += 1

def elapsed_time_graph(parameter, parameter_values, sets):
    with PdfPages('elapsed_time_graph_%s.pdf' % parameter) as pdf:
        os.chdir("data")
        dict = {}
        times = []
        for i in sets:
            os.chdir("set%d" % i)
            file = open("elapsed_time.txt")
            s = file.read()
            times.append(float(s))
            file.close()
            os.chdir("..")
        os.chdir("..")
        averages = []
        averages.append(sum(times[0:10])/10)
        averages.append(sum(times[10:20])/10)
        averages.append(sum(times[20:30])/10)
        # print(parameter_values)
        # print(averages)
        plt.plot(parameter_values,averages, "ro")
        plt.xlabel("%s values" % parameter)
        plt.ylabel("Elapsed Time")
        plt.title("Elapsed Time vs. %s" % parameter)
        pdf.savefig()
        plt.close()

def entropy_graph(parameter, parameter_values, sets):
    with PdfPages('entropy_graph_%s.pdf' % parameter) as pdf:
        os.chdir("data")
        dict = {}
        entropies = []
        for i in sets:
            os.chdir("set%d" % i)
            file = open("entropy.txt")
            s = file.read()
            entropies.append(float(s))
            file.close()
            os.chdir("..")
        os.chdir("..")
        averages = []
        averages.append(sum(entropies[0:10])/10)
        averages.append(sum(entropies[10:20])/10)
        averages.append(sum(entropies[20:30])/10)
        # print(parameter_values)
        # print(averages)
        plt.plot(parameter_values,averages, "ro")
        plt.xlabel("%s values" % parameter)
        plt.ylabel("Entropy")
        plt.title("Entropy vs. %s" % parameter)
        pdf.savefig()
        plt.close()

def position_overlap_graph(parameter, parameter_values, sets):
    with PdfPages('position_overlap_graph_%s.pdf' % parameter) as pdf:
        os.chdir("data")
        dict = {}
        pos_overlaps = []
        for i in sets:
            os.chdir("set%d" % i)
            file = open("overlap_pos.txt")
            s = file.read()
            pos_overlaps.append(float(s))
            file.close()
            os.chdir("..")
        os.chdir("..")
        averages = []
        averages.append(sum(pos_overlaps[0:10])/10)
        averages.append(sum(pos_overlaps[10:20])/10)
        averages.append(sum(pos_overlaps[20:30])/10)
        # print(parameter_values)
        # print(averages)
        plt.plot(parameter_values,averages, "ro")
        plt.xlabel("%s values" % parameter)
        plt.ylabel("Position Overlaps")
        plt.title("Position Overlap vs. %s" % parameter)
        pdf.savefig()
        plt.close()

def overlap_site_graph(parameter, parameter_values, sets):
    with PdfPages('overlap_site_graph_%s.pdf' % parameter) as pdf:
        os.chdir("data")
        dict = {}
        over_sites = []
        for i in sets:
            os.chdir("set%d" % i)
            file = open("overlap_sites.txt")
            s = file.read()
            over_sites.append(float(s))
            file.close()
            os.chdir("..")
        os.chdir("..")
        averages = []
        averages.append(sum(over_sites[0:10])/10)
        averages.append(sum(over_sites[10:20])/10)
        averages.append(sum(over_sites[20:30])/10)
        # print(parameter_values)
        # print(averages)
        plt.plot(parameter_values,averages, "ro")
        plt.xlabel("%s values" % parameter)
        plt.ylabel("Overlap Sites")
        plt.title("Overlap Site vs. %s" % parameter)
        pdf.savefig()
        plt.close()

def main():
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

        for i in range(1, 11):
            benchmark_helper(default_icpc, default_ml, default_sc)

        range_start = set_count
        for icpc in icpc_vals:
            for i in range(1, 11):
                benchmark_helper(icpc, default_ml, default_sc)

        sets = list(range(1, 11)) + list(range(range_start, set_count))
        elapsed_time_graph("ICPC", [default_icpc] + icpc_vals, sets)
        entropy_graph("ICPC", [default_icpc] + icpc_vals, sets)
        position_overlap_graph("ICPC", [default_icpc] + icpc_vals, sets)
        overlap_site_graph("ICPC", [default_icpc] + icpc_vals, sets)

        range_start = set_count
        for ml in motif_lengths:
            for i in range(1, 11):
                benchmark_helper(default_icpc, ml, default_sc)

        sets = list(range(1, 11)) + list(range(range_start, set_count))
        elapsed_time_graph("ML", [default_ml] + motif_lengths, sets)
        entropy_graph("ML", [default_ml] + motif_lengths, sets)
        position_overlap_graph("ML", [default_ml] + motif_lengths, sets)
        overlap_site_graph("ML", [default_ml] + motif_lengths, sets)

        range_start = set_count
        for sc in sequence_counts:
            for i in range(1, 11):
                benchmark_helper(default_icpc, default_ml, sc)

        sets = list(range(1, 11)) + list(range(range_start, set_count))
        elapsed_time_graph("SC", [default_sc] + sequence_counts, sets)
        entropy_graph("SC", [default_sc] + sequence_counts, sets)
        position_overlap_graph("SC", [default_sc] + sequence_counts, sets)
        overlap_site_graph("SC", [default_sc] + sequence_counts, sets)

if __name__ == "__main__":
    main()
