import math

import scipy.misc as spm
import scipy.special as sps

log = math.log
sqrt = math.sqrt


# A beta negative binomial distributed random variable
class BetaNegativeBinomial:
    def __init__(self, n, alpha, beta):
        self.n = n
        self.alpha = alpha
        self.beta = beta

    # Probability mass function
    def pmf(self, k):
        r, alpha, beta = self.n, self.alpha, self.beta

        return sps.gamma(r + k) / (
            math.factorial(k) * sps.gamma(r)) * sps.beta(
                alpha + r, beta + k) / sps.beta(alpha, beta)

    # Inverse cumulative density function that sums the true PDF
    def iCDF(self, p0):
        k = 0
        S = 0.0

        while S < p0:
            P = self.pmf(k)
            S += P
            k += 1

        return k - 1

    # Inverse cumulative density function that uses Stirling's approximation and
    # some term rewriting
    def iCDFapprox(self, p0):
        n, alpha, beta = self.n, self.alpha, self.beta
        k = 0
        S = 0.0
        t = 1

        ab = alpha + beta
        an = alpha + n

        bound = p0 * sqrt(ab / alpha / beta) * pow(
            ab, alpha * (log(alpha, ab) - 1) + beta * (log(beta, ab) - 1))
        upperbound = 4 * n

        while S < bound and k < upperbound:
            bi = beta + k
            anbi = alpha + n + beta + k
            P = t * sqrt(anbi / an / bi) * pow(
                anbi, an * (log(an, anbi) - 1) + bi * (log(bi, anbi) - 1))
            S += P
            k += 1
            t *= (n + k - 1) / (k)

        return k - 1
