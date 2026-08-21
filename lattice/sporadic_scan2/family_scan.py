#!/usr/bin/env python3
"""Zagier-style integrality search in the normalisation suggested by the two known
Sym^1 rows (Domb curve: (A,C)=(20,16); T curve: (A,C)=(24,4)):

    (n+1)^2 a_{n+1} = (A n^2 + (A/2) n + B) a_n - C (2n-1)^2 a_{n-1},   a_0=1, a_{-1}=0
    => a_1 = B.

This is NOT Zagier's class (his p_0 = c n^2 vanishes at n=0); it is the class the
weight-one CM rows live in.  Scan (A,B,C) and keep those with a_n in Z.
"""
from fractions import Fraction
import sys, json

NTEST = 60
def integral(A, B, C, nmax=NTEST):
    a0, a1 = Fraction(1), Fraction(B)
    if a1.denominator != 1: return None
    seq = [a0, a1]
    for n in range(1, nmax):
        num = Fraction(A*n*n) + Fraction(A, 2)*n + B
        val = (num*seq[n] - Fraction(C*(2*n-1)**2)*seq[n-1]) / Fraction((n+1)**2)
        if val.denominator != 1: return None
        seq.append(val)
    return [int(x) for x in seq]

out = []
import os
AL = int(os.environ.get("AL", "200")); CL = int(os.environ.get("CL", "120")); BL = int(os.environ.get("BL", "24"))
for A in range(-AL, AL+1, 2):              # A even so that A/2 in Z
    for B in range(-BL, BL+1):
        for C in range(-CL, CL+1):
            s = integral(A, B, C, 12)
            if s is None: continue
            s = integral(A, B, C, NTEST)
            if s is None: continue
            out.append({"A": A, "B": B, "C": C, "a": [str(x) for x in s[:10]]})
print(len(out), "integral triples")
json.dump(out, open("family_hits.json", "w"))
for h in out:
    print(h["A"], h["B"], h["C"], ",".join(h["a"][:7]))
