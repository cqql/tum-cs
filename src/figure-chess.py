#!/usr/bin/env python

import argparse

import data
import matplotlib.pyplot as pp
import numpy as np
import scipy
import sdpclstr


def solve(models):
    k = len(models)
    P = np.hstack([m.samples() for m in models])

    # P = P.T
    # np.random.shuffle(P)
    # P = P.T

    m, N = P.shape

    X_D = sdpclstr.sdp(P, k)

    return X_D


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--seed", default=10, type=int, help="Random seed")

    args = parser.parse_args()
    np.random.seed(args.seed)

    scenarios = [
        # 3 Normal Balls
        [data.BallCluster(10, [-2, -1], 1),
         data.BallCluster(15, [-1, 2.5], 1),
         data.BallCluster(15, [3, 1], 1.5)],
        # 3 Normal Balls
        [data.NormalCluster(10, [-2, -1], 1),
         data.NormalCluster(15, [-1, 2.5], 1),
         data.NormalCluster(15, [3, 1], 1.5)]
    ] # yapf: disable

    X_Ds = list(map(solve, scenarios))

    # Add an exact solution
    X_R = scipy.linalg.block_diag(
        np.ones((10, 10)) / 10, np.ones((15, 15)) / 5, np.ones((15, 15)) / 15)
    X_Ds.insert(0, X_R)

    logX_Ds = list(map(np.log10, X_Ds))

    fig, axes = pp.subplots(nrows=1, ncols=3)
    fig.set_size_inches(14, 4)

    for i, ax in enumerate(axes.flat):
        img = ax.matshow(
            logX_Ds[i], cmap=pp.get_cmap("Greys"), vmin=-15.0, vmax=0.0)

    fig.subplots_adjust(right=0.87)
    cbar = fig.add_axes([0.9, 0.12, 0.03, 0.76])
    fig.colorbar(img, cax=cbar)
    # fig.show()
    fig.savefig("../figures/chessboard.eps", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
