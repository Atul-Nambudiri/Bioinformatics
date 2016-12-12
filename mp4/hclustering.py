from math import sqrt
import numpy as np

def distance(a1, b1):
    min_dist = 10000000000
    for a in a1:
        for b in b1:
            min_dist = min(min_dist, sqrt((a[0] - b[0])**2 + (a[1] - b[1]) ** 2))
    return min_dist

def main():
    clusters = [
            [(1.0, 2.5)],
            [(1.0, 2.0)],
            [(3.0, 2.0)],
            [(3.0, 1.0)],
            [(3.5, 2.5)]]
    
    """clusters = [
            [(1.0, 2.5)],
            [(2.0, 3.0)],
            [(3.0, 1.0)],
            [(3.0, 2.0)],
            [(3.5, 2)]]
    """
    while(len(clusters) > 1):
        print("\n".join(str(cluster) for cluster in clusters))
        print()
        cluster_matrix = np.zeros((len(clusters), len(clusters)))
        for i in range(len(clusters)):
            for j in range(len(clusters)):
                if i == j:
                    cluster_matrix[i][j] = 100000000000000
                else:
                    cluster_matrix[i][j] = distance(clusters[i], clusters[j])
        pos = np.argmin(cluster_matrix)
        posi = pos / len(clusters)
        posj = pos % len(clusters)
        new_clusters = []
        for i in range(len(clusters)):
            if i == posi:
                new_clusters.append(clusters[posi] + clusters[posj])
            elif i != posj:
                new_clusters.append(clusters[i])
        clusters = new_clusters
    print("\n".join(str(cluster) for cluster in clusters))
    print("")

if __name__ == "__main__":
    main()
