from math import sqrt

def distance(a, b):
    return sqrt((a[0] - b[0])**2 + (a[1] - b[1]) ** 2)

def main():
    points = [(1.0, .5),
            (2.5, 3.0),
            (2.0, 1.0),
            (3.0, 2.0),
            (3.5, 2.0),
            (-0.5, 0.0),
            (-0.5, 1.0),
            (-1.0, 0.5),
            (1.0, -1.0),
            (0.5, 1.0)]
    k = 3
    cluster_centers = [(0.0, 0.0), (2.0, 3.0), (1.5, -1.0)]
    iteration = 0
    while(True):
        cluster_assignments = {0 : [], 1 : [], 2 : []}
        for point in points:
            cur_min = 10000000
            min_num = -1
            for i, cluster in enumerate(cluster_centers):
                dist = distance(point, cluster)
                if dist < cur_min:
                    min_num = i
                    cur_min = dist
            cluster_assignments[min_num].append(point)
        new_centers = []
        print("Step: %d" % iteration)
        for r, t in cluster_assignments.items():
            print("Cluster Center: (%3.4f; %3.4f)" % (cluster_centers[r][0], cluster_centers[r][1]))
            print("Points, Distance to Center\n" + "\n".join("(%3.4f; %3.4f), %3.4f" % (p[0], p[1], distance(p, cluster_centers[r])) for p in t))
            print("")
            top = sum(j[0] for j in t)
            bottom = sum(j[1] for j in t)
            new_centers.append((top/3.0, bottom/3.0))
        if cmp(cluster_centers, new_centers) == 0:
            return
        else:
            cluster_centers = new_centers
        iteration += 1
     




if __name__ == "__main__":
    main()
