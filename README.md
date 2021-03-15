# Module - netpy      version: 1.0
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
