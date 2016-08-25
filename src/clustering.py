#!/usr/bin/env python

import argparse

import data
import matplotlib.pyplot as pp
import numpy as np
import scipy
import sdpclstr


def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cluster",
        default=[data.BallCluster(15, [-2, -1], 1),
                 data.BallCluster(15, [3, 1], 1.5)],
        type=data.fromfile,
        nargs="*",
        help="Cluster specifications")
    parser.add_argument(
        "-e", "--epsilon", default=None, type=float, help="Epsilon parameter")
    parser.add_argument("--seed", default=10, type=int, help="Random seed")

    args = parser.parse_args()

    # Seed random state
    np.random.seed(args.seed)

    # Generate points
    models = args.cluster
    k = len(models)
    samples = [m.samples() for m in models]
    P = np.hstack(samples)
    m, N = P.shape

    # Save true cluster centers
    truecenters = [m.center() for m in models]

    # Save the k-means-optimal centroids
    optimalcentroids = [np.mean(S, axis=1) for S in samples]

    # Save true cluster affinities
    offset = 0
    trueclusters = []
    for S in samples:
        indices = set(range(S.shape[1]))
        indinP = set(map(lambda i: i + offset, indices))

        trueclusters.append(indinP)
        offset += S.shape[1]

    # Clustering
    PX_D, clusters, centers = sdpclstr.cluster(k, P, args.epsilon)

    ### Plotting ###

    # Classification error
    misclass = []
    tmp = clusters.copy()
    for C in trueclusters:
        # Match each cluster with the set of points that has the largest
        # intersection with it.
        match = max(tmp, key=lambda c: len(c.intersection(C)))

        misclass += list(match.difference(C))

        # Every cluster can only be matched once
        tmp.remove(match)

    error = len(misclass) / N
    msg = "Correctly classified {} / {}; Error = {:.3f}"
    print(msg.format(N - len(misclass), N, error))

    # Clustering error (distance to true centers)
    trueerror = 0
    optimalerror = 0
    tmp = set(range(centers.shape[1]))
    for i in range(k):
        tcenter = truecenters[i]
        optimal = optimalcentroids[i]
        match = min(tmp,
                    key=lambda c: scipy.linalg.norm(centers[:, c] - tcenter))

        trueerror += scipy.linalg.norm(centers[:, match] - tcenter)
        optimalerror += scipy.linalg.norm(centers[:, match] - optimal)

        tmp.remove(match)

    print("Error (True Centers) = {:.3f}".format(trueerror))
    print("Error (k-means-optimal centroids) = {:.3f}".format(optimalerror))

    if m != 2:
        return

    # Plot original and denoised data colored by cluster
    COLORS = ["c", "m", "b", "g", "y"]
    for i, cluster in enumerate(clusters):
        orig = P[:, list(cluster)]
        denoised = PX_D[:, list(cluster)]

        pp.plot(
            orig[0, :],
            orig[1, :],
            "o",
            color=COLORS[i],
            label="Cluster {}".format(i + 1))

        pp.plot(
            denoised[0, :],
            denoised[1, :],
            "d",
            ls="",
            ms=12,
            color=COLORS[i],
            label="Denoised Data {}".format(i + 1))

        # Plot k-means-optimal centroids
        centroid = optimalcentroids[min(len(models) - 1, i)]
        pp.plot(
            centroid[0],
            centroid[1],
            "h",
            ms=12,
            ls="",
            color=COLORS[i],
            label="Optimal Centroid {}".format(i + 1))

        # Plot true cluster centers
        center = truecenters[min(len(models) - 1, i)]
        pp.plot(
            center[0],
            center[1],
            "^",
            ms=15,
            ls="",
            color=COLORS[i],
            label="True Center {}".format(i + 1))

    # Mark misclassified points
    Pmisclass = P[:, misclass]
    pp.plot(
        Pmisclass[0, :], Pmisclass[1, :], "xr", ms=12, label="Misclassified")

    # Plot cluster center approximators
    pp.plot(
        centers[0, :], centers[1, :], "*y", label="Approximate centers", ms=12)

    pp.legend(loc="best", framealpha=0.5, numpoints=1)
    pp.show()


if __name__ == "__main__":
    main()
