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

        bound = p0 * sqrt(ab * an / alpha / beta)
        upperbound = 4 * n

        while S < bound and k < upperbound:
            bi = beta + k
            anbi = an + bi

            P = t * pow(an * ab / anbi / alpha, alpha) * pow(
                an / anbi,
                n) * pow(bi * ab / anbi / beta, beta) * pow(bi / anbi, k - 0.5)

            S += P
            k += 1
            t *= (n + k - 1) / (k)

        return k - 1
