"""Collect all optimisation results, recompute at high resolution, and emit the
final per-target table (including the irrationality measure kappa)."""
import json, glob, math, os, sys, warnings
import numpy as np
warnings.filterwarnings('ignore')
import haupt, outer, family, freeopt
from targets import TARGETS


def kappa(T, margin):
    """CDT ICM 6.2 measure: m = A/(B - logR * (2(1-g)/k - (1-g)^2/k^2)),
    g = 1/m.  Solve for kappa."""
    m, L = T['m'], T['L']
    g = 1.0 / m
    D = margin / m
    if margin <= 0:
        return None
    s = 1.0 - D / L
    if s < 0:
        return None
    den = 1.0 - math.sqrt(s)
    if den <= 0:
        return None
    return (1 - g) / den


def full_eval(T, u, Ns=(4096, 8192, 16384)):
    H = T['H']
    rows = []
    for N in Ns:
        uu = np.interp(np.arange(N) / N, np.arange(len(u)) / len(u), u, period=1.0) if len(u) != N else u
        r = outer.evaluate(H, uu, True)
        if r is None:
            continue
        rows.append((N, r))
    return rows


def report_row(T, V, lg):
    den = lg + T['L'] - T['tau']
    num = V + T['L']
    return dict(val=V, logdr=lg, den=den, num=num, cost=V - T['m'] * lg,
                margin=T['m'] * den - num, bound=(num / den if den > 1e-12 else None))
