import numpy as np
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

    def center(self):
        return self.p


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

    def center(self):
        return self.mu


def fromfile(path):
    with open(path, "r") as f:
        return eval(f.read(), {"__builtins__": None}, {"Ball": BallCluster,
                                                       "Normal":
                                                       NormalCluster})
