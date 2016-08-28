#!/usr/bin/env python

import matplotlib.pyplot as pp
import numpy as np
import scipy.stats
import sdpclstr


def main():
    P = np.vstack(
        (scipy.stats.multivariate_normal.rvs(mean=[-1, 2], cov=1, size=20),
         scipy.stats.multivariate_normal.rvs(mean=[-1, -1], cov=1, size=20),
         scipy.stats.multivariate_normal.rvs(mean=[2, 0], cov=1, size=20))).T
    Y = np.array([[-1, -1, 2], [2, -1, 0]])

    X_D = sdpclstr.sdp(P, 3)
    PX_D = P @ X_D

    # fig, ax = pp.subplots(nrows=1, ncols=2)
    # fig.set_size_inches(8, 4)
    fig1 = pp.figure(figsize=(4, 4))
    fig2 = pp.figure(figsize=(4, 4))

    # Set axis limits
    fig1.gca().set_xlim(-4, 4)
    fig2.gca().set_xlim(-4, 4)
    fig1.gca().set_ylim(-4, 4)
    fig2.gca().set_ylim(-4, 4)

    # Set ticks
    fig1.gca().set_xticks([-4, -2, 0, 2, 4])
    fig1.gca().set_yticks([-4, -2, 0, 2, 4])
    fig2.gca().set_xticks([-4, -2, 0, 2, 4])
    fig2.gca().set_yticks([-4, -2, 0, 2, 4])

    # Plot centers
    fig1.gca().plot(Y[0, :], Y[1, :], "xk", ms=15)
    fig2.gca().plot(Y[0, :], Y[1, :], "xk", ms=15)

    # Plot samples
    fig1.gca().plot(P[0, :], P[1, :], "ow")
    fig2.gca().plot(PX_D[0, :], PX_D[1, :], "ow")

    # fig.show()
    fig1.savefig("../figures/denoising-noisy.eps", dpi=160)
    fig2.savefig("../figures/denoising-denoised.eps", dpi=160)

if __name__ == "__main__":
    main()
