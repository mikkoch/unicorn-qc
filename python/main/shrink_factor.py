import numpy as np
import sympy as sym
from scipy.optimize import fsolve


def mp_cdf(x, rho, sigma=1.0):
    a = (1 - rho**0.5)**2
    b = (1 + rho**0.5)**2
    y = x / sigma
    int1 = np.arcsin((2 * y - a - b) / (b - a))
    int2 = np.arcsin(((a + b) * y - 2 * a * b) / ((b - a) * y)) / (a * b)**0.5
    return 0.5 + (((y - a) * (b - y))**0.5 + (a + b) * int1 / 2 - a * b * int2) / (2 * np.pi * rho)


def double_solve(f1, f2, x0, y0):
    func = lambda x: [f1(x[0], x[1]), f2(x[0], x[1])]
    return fsolve(func,[x0, y0])


def fit_median_iqr(eigs, m_guess, sigma_guess):
    med = np.median(eigs)
    n_samples = len(eigs)
    q75, q25 = np.percentile(eigs, [75 ,25])
    median = lambda m, s: mp_cdf(med, n_samples/m, s) - 0.5
    iqr = lambda m, s: mp_cdf(q75, n_samples/m, s) - mp_cdf(q25, n_samples/m, s) - 0.5
    return double_solve(median, iqr, m_guess, sigma_guess)


def shrinkage_scores(eig, n_samples, n_variants):
    x = sym.Symbol('x')
    gamma = n_samples / n_variants
    if eig > (1 + gamma**0.5)**2:
        lmbd = sym.solve(x + gamma * x / (x - 1) - eig, x)[1]
    else: lmbd = 1
    return (lmbd - 1) / lmbd