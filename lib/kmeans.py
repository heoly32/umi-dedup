from itertools import imap, combinations
import collections, operator, random

def hamming(umi1, umi2):
        assert len(umi1) == len(umi2)
        ne = operator.ne
        return sum(imap(ne, umi1, umi2))

class DistanceMatrix:
    '''Create object that store distance matrices'''

    def __init__(self, graph, counts):
        self.distances = {(umi1, umi2): min(counts[umi1], counts[umi2])/max(counts[umi1], counts[umi2]) * hamming(umi1.encode('utf-8'),
                                                                                                                  umi2.encode('utf-8')) for umi1, umi2 in combinations(graph.keys(), 2)}

    def get_distance(self, umi1, umi2):
        if umi1 == umi2:
            return 0
        else:
            try:
                return self.distances[(umi1, umi2)]
            except KeyError:
                return self.distances[(umi2, umi1)]

def cluster_points(distance_matrix, umis, centroids):
    clusters  = {}
    for umi in umis:
        distances = [distance_matrix.get_distance(center, umi) for center in centroids]
        closest_centroid = distances.index(min(distances))
        try:
            clusters[centroids[closest_centroid]].add(umi)
        except KeyError:
            clusters[centroids[closest_centroid]] = set([umi])
    return clusters

def reevaluate_centers(centroids, clusters, counts):
    new_centroids = []
    for cluster in clusters.values():
        cluster = list(cluster)
        counts_restricted = [counts[umi] for umi in cluster]
        most_abundant = counts_restricted.index(max(counts_restricted))
        new_centroids.append(cluster[most_abundant])
    return new_centroids

def has_converged(centroids, old_centroids):
    return set(centroids) == set(old_centroids)

def find_clusters(graph, counts, k):
    # Initialize to k random centers
    umis = graph.keys()
    old_centroids = random.sample(graph, k)
    centroids = random.sample(graph, k)
    clusters = {"graph": set(umis)}

    # print centroids, old_centroids

    distance_matrix = DistanceMatrix(graph, counts)

    while not has_converged(centroids, old_centroids):
        old_centroids = centroids
        # Assign all points in X to clusters
        clusters = cluster_points(distance_matrix, umis, centroids)
        # Reevaluate centers
        centroids = reevaluate_centers(old_centroids, clusters, counts)
        # print centroids, old_centroids

    return clusters
