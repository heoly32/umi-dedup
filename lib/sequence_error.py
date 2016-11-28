from itertools import imap
import collections, operator

def breadth_first_search(node, adj_list):
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))

    while len(queue) > 0:
        node = (list(queue))[0]
        found.update(adj_list[node])
        queue.update(adj_list[node])
        searched.update((node,))
        queue.difference_update(searched)

    return found

def remove_umis(adj_list, cluster, nodes):
    '''removes the specified nodes from the cluster and returns
    the remaining nodes '''

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove

def hamming(umi1, umi2):
    assert len(umi1) == len(umi2)
    ne = operator.ne
    return sum(imap(ne, umi1, umi2))

class ClusterAndReducer:
    '''A functor that clusters a bundle of reads,
    indentifies the parent UMIs and returns the selected reads, umis and counts
    Methods:
      ** get_adj_list ** - returns the edges connecting the UMIs
      ** connected_components ** - returns clusters of connected components
                                   using the edges in the adjacency list
      ** get_best ** - returns the parent UMI(s) in the connected_components
      ** reduce_clusters ** - loops through the connected components in a
                              cluster and returns the unique reads. Optionally
                              returns lists of umis and counts per umi also
    '''

    def get_best(self, cluster, counts):
        ''' return the UMI with the highest counts'''
        if len(cluster) == 1:
            return list(cluster)[0]
        else:
            sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                                  reverse=True)
            return sorted_nodes[0]

    def get_adj_list(self, umis, counts, threshold = 1):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (2 * second umi counts)-1'''

        return {umi: [umi2 for umi2 in umis if
                      hamming(umi.encode('utf-8'),
                                    umi2.encode('utf-8')) == threshold and
                      counts[umi] >= (counts[umi2] * 2) - 1] for umi in umis}

    def get_connected_components(self, umis, graph, counts):
        ''' find the connected UMIs within an adjacency dictionary'''

        found = list()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)

        return components

    def reduce_clusters(self, bundle, clusters,
                                  adj_list, counts):
        ''' collapse clusters down to the UMI which accounts for the cluster
        using the adjacency dictionary and return the list of final UMIs'''

        reads = {}

        for cluster in clusters:
            parent_umi = get_best(cluster, counts)
            reads[parent_umi] = [read for umi in cluster for read in bundle[umi]]

        return reads

    def __call__(self, bundle, threshold = 1):

        umis = bundle.keys()

        len_umis = [len(x) for x in umis]
        assert max(len_umis) == min(len_umis), (
            "not all umis are the same length(!):  %d - %d" % (
                min(len_umis), max(len_umis)))

        counts = {umi: len(bundle[umi]) for umi in umis}

        adj_list = self.get_adj_list(umis, counts, threshold)

        clusters = self.get_connected_components(umis, adj_list, counts)

        reads = self.reduce_clusters(bundle, clusters, adj_list, counts)

        return reads
