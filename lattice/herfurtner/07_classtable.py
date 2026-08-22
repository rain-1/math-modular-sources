#!/usr/bin/env python3
"""Normalisation classes <-> Herfurtner configurations.

A second-order Apery-like row with four singular points has, at 0, a cusp (MUM);
at t_1,t_2 (Galois-conjugate unless rational) equal PSL_2(Z) local type; at oo a
third type.  Group the rigid Herfurtner configurations by the pair
(type at t_1,t_2 ; type at oo) and list the cross-ratio invariant I = A^2/D.
"""
import sys
from sympy import simplify, N as num, Rational as Q
sys.path.insert(0,'.')
from importlib import import_module
H = import_module('04_herfurtner')

T = H.table()
groups = {}
for r in T:
    key = (r['rho'], r['dinf'])
    groups.setdefault(key, {})
    k2 = (r['row'], r['at0'], r['atinf'], r['att'])
    groups[key][k2] = r

NAME = {'c':'cusp (I_n / I_n^*)', 2:'order 2 (III/III^*)', 3:'order 3 (II,IV,IV^*,II^*)'}
QCLASS = {('c','c'):'(1;0,0) rho=0 delta=0  [Zagier]  and (1;-1,1),(2;-1,1),(1;1,1),(2;1,3)',
          ('c',2):'(4;-1,1) rho=0 delta=1/2 ; (4;3,5) rho=1 delta=1/2',
          ('c',3):'(3;-1,1) rho=0 delta=2/3 ; (6;-1,1) rho=0 delta=1/3 ; (3;2,4),(6;5,7)',
          (2,'c'):'(2;1,1) rho=1/2 delta=0  [root rows] ; (1;0,1) rho=1/2 delta=1',
          (2,2):'(4;1,3) rho=1/2 delta=1/2',
          (2,3):'(3;1,2) rho=1/2 delta=1/3 ; (6;1,5) rho=1/2 delta=2/3',
          (3,'c'):'(3;1,1) rho=1/3 ; (3;2,2) rho=2/3',
          (3,2):'(12;1,7) rho=1/3 delta=1/2 ; (12;5,11) rho=2/3 delta=1/2',
          (3,3):'(3;0,2),(3;1,3) delta=2/3 ; (6;1,3),(6;3,5) delta=1/3'}

for key in sorted(groups, key=lambda k:(str(k[0]),str(k[1]))):
    tt, ti = key
    print("="*100)
    print("t_1,t_2 : %-28s   infinity : %-28s" % (NAME[tt], NAME[ti]))
    print("recurrence classes: %s" % QCLASS.get(key,'?'))
    seen=set()
    for k2,r in sorted(groups[key].items(), key=lambda kv: kv[0][0]):
        Iv = r['I']
        try: f = complex(num(Iv))
        except Exception: f = None
        tag = "%.6f"%f.real if (f is not None and abs(f.imag)<1e-20) else str(f)
        print("   #%-3d %-22s  degJ=%-3d  0:%-5s oo:%-5s  t:%s,%s   I = %-28s ~ %s"
              % (r['row']," ".join(r['fib']),r['degJ'],r['at0'],r['atinf'],
                 r['att'][0],r['att'][1], str(Iv), tag))
