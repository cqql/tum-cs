#!/usr/bin/env python

import math

import numpy as np
import matplotlib.pyplot as pp


def setupfig():
    pp.ylim((0, 1))
    ax = pp.gca()
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.xaxis.set_ticks([0, 1])
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks([0, 0.5, 1])
    ax.yaxis.set_ticks_position("left")

    pp.gcf().set_size_inches((4, 4))


def main():
    X = np.linspace(0, 1, 200)

    Y = 0.5 * np.ones((200, ))
    pp.plot(X, Y, "k--", lw=2)
    setupfig()
    pp.savefig("constant.eps")
    pp.clf()

    Y = list(map(lambda x: x, X))
    pp.plot(X, Y, "k--", lw=2)
    setupfig()
    pp.savefig("linear.eps")
    pp.clf()

    def f(x):
        if 0.25 <= x <= 0.75:
            return math.cos(4 * (x - 0.25) * math.pi - math.pi / 2) / 4 + 0.5
        else:
            return x

    Y = list(map(f, X))
    pp.plot(X, Y, "k--", lw=2)
    setupfig()
    ax = pp.gca()
    ax.xaxis.set_ticks([0, 3/8, 5/8, 1])
    ax.yaxis.set_ticks([0, 0.75, 0.25, 1])
    pp.grid()
    pp.savefig("higher.eps")
    pp.clf()


if __name__ == "__main__":
    main()
