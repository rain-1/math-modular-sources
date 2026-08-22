#!/usr/bin/env python3
import sys, json
from sympy import Rational as Q, simplify
sys.path.insert(0,'.')
from importlib import import_module
H = import_module('04_herfurtner')

T = H.table()

def show(name,M,j1,j2,A,B,C):
    Iv,(tt,ti),hits = H.match_row(M,j1,j2,A,C,T)
    D=C*M*M
    s = "%-24s class(%d;%d,%d) A=%-5d B=%-5d D=%-7d  I=%-12s types(t,oo)=(%s,%s)  " % (
        name,M,j1,j2,A,B,D,str(Iv),tt,ti)
    if hits:
        s += " | ".join("#%d %s [0:%s, oo:%s, t:%s,%s] degJ=%d"%(r['row']," ".join(r['fib']),
              r['at0'],r['atinf'],r['att'][0],r['att'][1],r['degJ']) for r in hits)
    else:
        s += "NO rigid Herfurtner configuration"
    print(s)

ROWS=[
 ("Zagier A (7,2,-8)",1,0,0,7,2,-8),
 ("Zagier B (9,3,27)",1,0,0,9,3,27),
 ("Zagier C (10,3,9)",1,0,0,10,3,9),
 ("Zagier D (11,3,-1)",1,0,0,11,3,-1),
 ("Zagier E (12,4,32)",1,0,0,12,4,32),
 ("Zagier F (17,6,72)",1,0,0,17,6,72),
 ("root Apery",2,1,1,136,10,4),
 ("root T",2,1,1,24,2,4),
 ("root Domb",2,1,1,20,2,16),
 ("root AZ(9,3,-27)",2,1,1,72,6,-108),
 ("root AZ(11,5,125)",2,1,1,88,10,500),
 ("root AZ(7,3,81)",2,1,1,56,6,324),
 ("root Cooper s7",3,1,2,26,2,-3),
]
for r in ROWS: show(*r)
