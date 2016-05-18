#!/usr/bin/env python

import argparse
import collections
import operator as op

import numpy as np

import cvxopt as cvx
import matplotlib.pyplot as pp
import networkx as nx
import picos as pic
import scipy
import scipy.spatial
import scipy.stats


# Stochastic ball model cluster
class BallCluster:
    def __init__(self, n, p, r):
        self.n = int(n)
        self.p = np.array(p)
        self.r = float(r)

    # Sample points in the r-ball with rejection sampling
    def samples(self):
        m = len(self.p)
        P = np.zeros((m, self.n))
        for i in range(self.n):
            p = scipy.stats.uniform.rvs(loc=-self.r, scale=2 * self.r, size=m)

            while scipy.linalg.norm(p) > self.r:
                p = scipy.stats.uniform.rvs(size=m)

            P[:, i] = p + self.p

        return P


# Normally distributed cluster
class NormalCluster:
    def __init__(self, n, mu, sigma):
        self.n = int(n)
        self.mu = np.array(mu)
        self.sigma = np.array(sigma)

    def samples(self):
        return scipy.stats.multivariate_normal.rvs(mean=self.mu,
                                                   cov=self.sigma,
                                                   size=self.n).T


def parsecluster(string):
    return eval(string, {"__builtins__": None}, {"Ball": BallCluster,
                                                 "Normal": NormalCluster})

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "--cluster",
    default=[BallCluster(15, [-2, -1], 1), BallCluster(15, [3, 1], 1.5)],
    type=parsecluster,
    nargs="*",
    help="Cluster specifications")
parser.add_argument("-e",
                    "--epsilon",
                    default=None,
                    type=float,
                    help="Epsilon parameter")
parser.add_argument("--seed", default=10, type=int, help="Random seed")

args = parser.parse_args()

### Setup ###

# Seed random state
np.random.seed(args.seed)

# Generate points
clusters = args.cluster
k = len(clusters)
P = np.hstack([c.samples() for c in clusters])
N = P.shape[1]

### Clustering ###

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

# Denoised data
PX_D = P @X_D

# Extract approximate cluster affinities
nodes = set(range(N))
clusters = []
while len(nodes) > 0:
    node = nodes.pop()
    row = X_D[node, :]

    cluster = [node]
    for i in range(N):
        if i in nodes and row[i] > 1 / N:
            cluster.append(i)
            nodes.remove(i)

    clusters.append(cluster)

# Rounding
if args.epsilon:
    epsilon = args.epsilon
else:
    # Compute cluster means from denoised data
    means = [np.mean(PX_D[:, c], axis=1) for c in clusters]

    # Compute minimum distance as approximation to d_min
    d_min = min(scipy.spatial.distance.pdist(means))

    # Use something slightly smaller than the upper bound d_min / 8
    epsilon = d_min / 10

    print("Using epsilon = {}".format(epsilon))

# Construct "adjacency" graph
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
        print("Epsilon {} is too large. Graph is empty in round {}.".format(
            epsilon, i + 1))
        exit()

    center = sorted(G.degree_iter(), key=op.itemgetter(1), reverse=True)[0][0]

    for j in G.nodes():
        if D[center, j] <= 4 * epsilon:
            G.remove_node(j)

    centers.append(P[:, center])

### Plotting ###

# Plot original and denoised data colored by cluster
COLORS = ["r", "g", "b"]
for i, cluster in enumerate(clusters):
    orig = P[:, cluster]
    denoised = PX_D[:, cluster]

    pp.plot(orig[0, :],
            orig[1, :],
            "o",
            color=COLORS[i],
            label="Cluster {}".format(i + 1))

    pp.plot(denoised[0, :],
            denoised[1, :],
            "x",
            ls="",
            color=COLORS[i],
            label="Denoised Data {}".format(i + 1))

# Plot cluster centers
for i, center in enumerate(centers):
    pp.plot([center[0]],
            [center[1]],
            "*",
            label="Center {}".format(i + 1),
            ms=15)

pp.legend(loc="best", framealpha=0.5, numpoints=1)
pp.show()
