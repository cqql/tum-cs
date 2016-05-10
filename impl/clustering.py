#!/usr/bin/env python

import argparse
import operator as op

import numpy as np

import cvxopt as cvx
import matplotlib.pyplot as pp
import networkx as nx
import picos as pic
import scipy.spatial

# Generate points
n = 15
k = 2
d_min = np.sqrt(18)
P = []
P.append(np.random.normal(-1, 1, (2, n)))
P.append(np.random.normal(3, 2, (2, n)))
P = np.hstack(P)
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

# Plot original data colored by cluster
COLORS = ["r", "g", "b"]
nodes = set(range(N))
clusternr = 0
while len(nodes) > 0:
    node = nodes.pop()
    row = X_D[node, :]

    cluster = [node]
    for i in range(N):
        if i in nodes and row[i] > 1 / N:
            cluster.append(i)
            nodes.remove(i)

    points = P[:, cluster]
    pp.plot(points[0, :],
            points[1, :],
            "o",
            color=COLORS[clusternr],
            label="Cluster {}".format(clusternr + 1))

    clusternr += 1

# Plot denoised data
PX_D = P @X_D
pp.plot(PX_D[0, :], PX_D[1, :], "x", ls="", label="Denoised Data")

# Rounding
epsilon = d_min / 9

# Construct "adjacency" graph
G = nx.Graph()
G.add_nodes_from(range(N))
for i in range(N):
    for j in range(i):
        if D[i, j] <= 2 * epsilon:
            G.add_edge(i, j)

# Extract nodes of maximum degree as cluster centers
for i in range(k):
    center = sorted(G.degree_iter(), key=op.itemgetter(1), reverse=True)[0][0]

    for j in G.nodes():
        if D[center, j] <= 4 * epsilon:
            G.remove_node(j)

    point = P[:, center]
    pp.plot([point[0]], [point[1]], "*", label="Center {}".format(i))

pp.legend(loc="best", framealpha=0.5)
pp.show()
