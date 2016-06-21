#!/usr/bin/env python

import argparse

import matplotlib.pyplot as pp
import numpy as np
from scipy.interpolate import spline


def plot(x, u):
    # A finer grid
    n = 300
    X = np.linspace(x[0], x[-1], n)

    # True function
    y = np.ones(X.shape)
    y[X > 0.5] = 0
    pp.plot(X, y, "k", label="exact")

    interp = spline(x, u, X, order=3)
    pp.plot(X, interp, "k--", lw=2, label="numerical")

    pp.legend(loc="upper right", frameon=False)

    ax = pp.gca()
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")


def parseoctavefile(path):
    with open(path, "r") as f:
        for l in f:
            if not l.startswith("#"):
                return np.fromstring(l, sep=" ")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("xfile", help="Path to x coordinates")
    parser.add_argument("ufile", help="Path to u values")
    parser.add_argument("out", help="Path to output file")

    args = parser.parse_args()

    x = parseoctavefile(args.xfile)
    u = parseoctavefile(args.ufile)

    plot(x, u)

    fig = pp.gcf()
    fig.set_size_inches(5, 5)
    fig.savefig(args.out, frameon=False, dpi=100)


if __name__ == "__main__":
    main()
