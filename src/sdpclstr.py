import operator as op

import cvxopt as cvx
import networkx as nx
import numpy as np
import picos as pic
import scipy
import scipy.spatial


def sdp(P, k):
    """Solve the SDP relaxation.

    Parameters
    ----------
    P : np.array
        m*N matrix of N m-dimensional points stacked columnwise
    k : int
        Number of clusters

    Returns
    -------
    The SDP minimizer S_D
    """
    N = P.shape[1]

    # Compute pairwise distances
    D = scipy.spatial.distance.pdist(P.T)
    D = scipy.spatial.distance.squareform(D)

    # Create the semi-definite program
    sdp = pic.Problem()
    Dparam = pic.new_param("D", D)
    X = sdp.add_variable("X", (N, N), vtype="symmetric")
    sdp.add_constraint(pic.trace(X) == k)
    sdp.add_constraint(X * "|1|({},1)".format(N) == 1)
    sdp.add_constraint(X > 0)
    sdp.add_constraint(X >> 0)
    sdp.minimize(Dparam | X)

    # Extract the SDP solution
    X_D = np.array(X.value)

    return X_D


def cluster(k, P, epsilon=None):
    """Run the SDP clustering algorithm.

    Parameters
    ----------
    k : int
        Number of clusters
    P : np.array
        m*N matrix of N m-dimensional points stacked columnwise
    epsilon : float, optional
        Radius of the ball used for medoid determination in the end

        If no value is given, it will be estimated from the denoised data.

    Returns
    -------
    (PX_D, clusters, centers)
    PX_D : np.array
        Denoised points
    clusters : list(list(int))
        Lists of indices of points that probably belong to the same cluster
    centers : np.array
        m*k matrix of k cluster-medoids
    """
    N = P.shape[1]
    X_D = sdp(P, k)

    # Denoised data
    PX_D = P @ X_D # yapf: disable

    # Extract approximate cluster affinities
    nodes = set(range(N))
    clusters = []
    while len(nodes) > 0:
        node = nodes.pop()
        row = X_D[node, :]

        cluster = set([node])
        for i in range(N):
            if i in nodes and row[i] > 1 / N:
                cluster.add(i)
                nodes.remove(i)

        clusters.append(cluster)

    # Rounding
    if not epsilon:
        # Compute cluster means from denoised data
        means = [np.mean(PX_D[:, list(c)], axis=1) for c in clusters]

        # Compute minimum distance as approximation to d_min
        d_min = min(scipy.spatial.distance.pdist(means))

        # Use something slightly smaller than the upper bound d_min / 8
        epsilon = d_min / 10

        print("Using epsilon = {:.3f}".format(epsilon))

    # Compute pairwise distances of denoised data
    D = scipy.spatial.distance.pdist(PX_D.T)
    D = scipy.spatial.distance.squareform(D)

    # Construct "adjacency" graph from denoised distances
    G = nx.Graph()
    G.add_nodes_from(range(N))
    for i in range(N):
        for j in range(i):
            if D[i, j] <= 2 * epsilon:
                G.add_edge(i, j)

    # Extract nodes of maximum degree as cluster centers
    centers = []
    for i in range(k):
        if G.size() == 0:
            print("Epsilon {} is too large. Graph is empty in round {}.".
                  format(epsilon, i + 1))
            exit()

        center = sorted(
            G.degree_iter(), key=op.itemgetter(1), reverse=True)[0][0]

        for j in G.nodes():
            if D[center, j] <= 4 * epsilon:
                G.remove_node(j)

        centers.append(center)

    centers = PX_D[:, centers]

    return (PX_D, clusters, centers)
