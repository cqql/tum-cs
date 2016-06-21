#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pp


def main():
    X = np.linspace(0, 1, 200)
    Y = np.linspace(-1, 1, 200)
    pp.plot(X, Y, "k--")

    ax = pp.gca()
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.xaxis.set_ticks([0, 0.5, 1])
    ax.xaxis.set_ticklabels(["$x_i$", r"$\frac{x_{i + 1} - x_{i}}{2}$",
                             "$x_{i + 1}$"])
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks([-1, 0, 1])
    ax.yaxis.set_ticks_position("left")
    # ax.grid(linestyle="dotted")

    pp.gcf().set_size_inches(2, 3.5)

    pp.savefig("xi.eps", dpi=100)


if __name__ == "__main__":
    main()
