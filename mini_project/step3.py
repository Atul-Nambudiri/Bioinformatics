import os

SET_COUNT = 70

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

def main():
    set_string = "data/set"

    motif_file = "motif.txt"
    predicted_motif_file = "predictedmotif.txt"

    ml_file = "motiflength.txt"

    sites_file = "sites.txt"
    predicted_sites_file = "predictedsites.txt"

    for i in range(1, 71):
        path = set_string + str(i)

        ml = parse_ml(get_path(path, ml_file))
        sites = parse_sites(get_path(path, sites_file))
        p_sites = parse_sites(get_path(path, predicted_sites_file))

        overlap_sites = compute_overlapping_sites(p_sites, sites)
        overlap_pos = compute_overlapping_pos(p_sites, sites, ml)

if __name__ == "__main__":
    main()
