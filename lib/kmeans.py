import numpy as np

def compute_dist_matrix(graph, counts):

    distance_matrix = {}
    for umi1 in graph.keys():
        for umi2 in graph[umi1]:
            distance_matrix[(umi1, umi2)] = 1 - min(counts[umi1], counts[umi2])/max(counts[umi1], counts[umi2])

    return distance_matrix

def cluster_points(distance_matrix, umis, centroids):
    clusters  = {}
    for umi in umis:
        closest_centroid = centroids.index(min([distance_matrix[(center, umi)] for center in centroids]))
        try:
            clusters[centroid[closest_centroid]].add(x)
        except KeyError:
            clusters[centroid[closest_centroid]] = set([x])
    return clusters

def reevaluate_centers(centroids, clusters):
    new_centroids = []
    for key in sorted(clusters.keys()):
        new_centroids.append(np.mean(clusters[k], axis = 0))
    return new_centroids

def has_converged(centroids, old_centroids):
    return (set(centroids) == set(old_centroids))

def find_clusters(graph, k):
    # Initialize to k random centers
    old_centroids = random.sample(graph, k)
    centroids = random.sample(graph, k)

    distance_matrix = compute_dist_matrix(graph)

    while not has_converged(centroids, old_centroids):
        old_centroids = centroids
        # Assign all points in X to clusters
        clusters = cluster_points(distance_matrix, centroids)
        # Reevaluate centers
        centroids = reevaluate_centers(old_centroids, clusters)
    return clusters
