# Fluctuation analysis: simulate cultures and estimate mutations per culture (m) three ways.
import math
import numpy as np


def grow_culture(n0, generations, mu, rng):
    """Synchronous doubling; each division makes a mutant with probability mu, and mutants grow too."""
    wild, mutants = n0, 0
    for _ in range(generations):
        new = rng.binomial(wild, mu)
        mutants = 2 * mutants + new
        wild = 2 * wild - new
    return mutants


def induced_culture(final_cells, rate, rng):
    # The alternative hypothesis: resistance is induced on the plate, so counts are Poisson.
    return rng.poisson(final_cells * rate)


def m_p0(counts):
    p0 = np.mean(np.asarray(counts) == 0)
    return -math.log(p0) if 0 < p0 < 1 else float("nan")


def m_lea_coulson(counts):
    r = float(np.median(counts))
    lo, hi = 1e-4, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if r / mid - math.log(mid) - 1.24 > 0:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def ld_pmf(m, r_max):
    """Ma-Sandri-Sarkar recursion for the Luria-Delbruck distribution, P(0..r_max)."""
    p = np.zeros(r_max + 1)
    p[0] = math.exp(-m)
    for r in range(1, r_max + 1):
        i = np.arange(r)
        p[r] = m / r * np.sum(p[:r] / (r - i + 1))
    return p


def m_mss_mle(counts, lo=0.01, hi=50.0):
    counts = np.asarray(counts)
    r_max = int(counts.max())

    def nll(m):
        return -np.sum(np.log(np.maximum(ld_pmf(m, r_max)[counts], 1e-300)))

    a, b = math.log(lo), math.log(hi)
    g = (math.sqrt(5) - 1) / 2
    for _ in range(60):
        c, d = b - g * (b - a), a + g * (b - a)
        if nll(math.exp(c)) < nll(math.exp(d)):
            b = d
        else:
            a = c
    return math.exp((a + b) / 2)
