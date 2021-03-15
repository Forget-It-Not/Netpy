'''
Module - netpy      version: 1.0
--------------------------------

Module for the manipulation of complex networks and calculation of several
measurements

Measures included:
    - Degrees and degree distributions
    - Nth moment of degree distributions
    - Average nearest neighbor degrees and degree distributions
    - Clustering coefficients
    - Geodesic matrix
    - Characteristic path length and efficiency

Networks are handled in a memory-efficient edgelist format, that is, as a list
of (ni, nj) tuples where ni and nj are nodes.
'''

# Auxiliar module import
import pandas as pd
import numpy as np

# Sample network initialization
def main():
    mynet = [(10,7), (7,3), (3,4), (10,4)]

if __name__ == '__main__':
    main()

def nodes(net):
    '''
    Returns list of nodes in network.
    '''
    nodes = []
    for edge in net:
        if edge[0] not in nodes:
            nodes.append(edge[0])
        if edge[1] not in nodes:
            nodes.append(edge[1])
    return nodes

def degree(net, node):
    '''
    Returns degree of node in network.
    '''
    deg = 0
    for edge in net:
        if node in edge:
            deg += 1
    return deg

def induced(net, subnodes):
    '''
    Returns induced subgraph of a list of nodes (subnodes) in network.

    Induced subgraph is defined as the subgraph of the network such that all
    links between pairs of nodes present in subnodes are also included in the
    subgraph.
    '''
    subnet = []
    for edge in net:
        if edge[0] in subnodes and edge[1] in subnodes:
            subnet.append(edge)
    return subnet

def neighbors(net, node):
    '''
    Returns list of nodes adjacent to a node.
    '''
    neighbors = []
    for edge in net:
        if node in edge:
            neighbor = edge[0] if edge[0] != node else edge[1]
            neighbors.append(neighbor)
    return neighbors

def degree_distribution(net, normalize=True):
    '''
    Returns degree distribution of network.

    Keyword arguments:
        - normalize=(True/False): whether to normalize frequencies
    '''
    degree_dist = []
    degrees = []
    nods = nodes(net)
    for node in nods:
        degrees.append(degree(net, node))
    for k in range(max(degrees)+1):
        degree_dist.append(degrees.count(k))
    if normalize:
        degree_dist = list(np.asarray(degree_dist)/len(nods))
    return degree_dist

def n_moment(degree_dist, n):
    '''
    Returns nth moment of a degree distribution.
    '''
    n_moment = 0
    for k in range(len(degree_dist)):
        n_moment += k**n*degree_dist[k]
    return n_moment

def avg_nn_degree(net, node):
    '''
    Returns average nearest neighbor degree of a node.
    '''
    neighs = neighbors(net, node)
    degree_sum = 0
    for neighbor in neighs:
        degree_sum += degree(net, neighbor)
    return (degree_sum)/(degree(net, node))

def avg_nn_degree_dist(net, normalize=True):
    '''
    Returns average nearest neighbor degree distribution of network.
    '''
    avg_nn_degree_dist = []
    values = []
    nods = nodes(net)
    for node in nods:
        values.append((degree(net, node), avg_nn_degree(net, node)))
    for k in range(max(values)[0]+1):
        avg = 0
        count = 0
        for val in values:
            if val[0] == k:
                avg += val[1]
                count += 1
        if normalize and count != 0:
            avg /= count
        avg_nn_degree_dist.append(avg)
    return avg_nn_degree_dist

def ci(net, node):
    '''
    Returns local clustering coefficient of a node in network.
    '''
    deg = degree(net, node)
    if deg in [0,1]:
        return 0
    neighs = neighbors(net, node)
    Gi = induced(net, neighs)
    edges = len(Gi)
    return ((2*edges) / (deg*(deg-1)))

def clustering_coefficient(net):
    '''
    Returns global (average) clustering coefficient of network.
    '''
    summ = 0
    nods = nodes(net)
    for node in nods:
        summ += ci(net, node)
    return summ/len(nods)

def geodesic(net):
    '''
    Returns geodesic matrix of network.

    Geodesic matrix is a NxN matrix where N is the nº of nodes in which the
    value at position ni,nj is the shortest distance (nº edges) between nodes
    ni and nj.
    '''
    nods = nodes(net)
    N = len(nods)
    t = {node:i for node, i in zip(nods, range(N))}
    geodesic_matrix = np.ones((N,N), dtype=int)*(-1)
    for node in nods:
        step = 1
        geodesic_matrix[t[node],t[node]] = 0
        sources = {node}
        while sources != set():
            new_sources = set()
            for source in sources:
                neighs = neighbors(net, source)
                for neighbor in neighs:
                    if geodesic_matrix[t[node], t[neighbor]] == -1:
                        geodesic_matrix[t[node], t[neighbor]] = step
                        new_sources = new_sources.union({neighbor})
            sources = new_sources
            step += 1
    return geodesic_matrix, t

def charact_path_length(geo):
    '''
    Returns characteristic path length of a network's geodesic matrix.
    '''
    N = len(geo)
    if len(geo[geo==-1]) != 0:
        raise Exception('Geodesic matrix contains disconnected components. ' +
                        'Isolate each component or measure efficiency instead.')
    summ = sum(sum(geo))
    return summ/(2*N*(N-1))

def efficiency(geo):
    '''
    Returns efficiency of a network's geodesic matrix.
    '''
    summ = 0
    N = len(geo)
    for i in range(len(geo)):
        for j in range(len(geo)):
            if i!=j and geo[i,j] != -1:
                summ += 1/geo[i,j]
    return summ/(2*N*(N-1))
