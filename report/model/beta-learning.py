#!/usr/bin/env python

import matplotlib.pyplot as pp
import numpy as np
import scipy.stats as sps

def plot(a, b, out):
    X = np.arange(0, 1, 0.01)
    Y = [sps.beta.pdf(x, a, b) for x in X]
    pp.plot(X, Y)
    pp.savefig(out)
    pp.clf()

def main():
    plot(1, 1, "beta-learning-1-1.eps")
    plot(2, 1, "beta-learning-2-1.eps")
    plot(3, 1, "beta-learning-3-1.eps")
    plot(3, 2, "beta-learning-3-2.eps")

if __name__ == "__main__":
    main()
