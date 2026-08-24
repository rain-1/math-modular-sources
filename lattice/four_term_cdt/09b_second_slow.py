#!/usr/bin/env python3
"""09b_second_slow.py -- 09_second.py for the two rows whose fold/second-singularity
ratio is close to 1 (R1: 0.828, R4: 0.933), where lim b^(k)_n/a_n needs many more
terms.  Same question: is xi^(2) (the fold constant of the solution of L y = t^2)
in Q + Q zeta(2) + Q G (resp. + Q L(2,chi_-3))?"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location("s9", os.path.join(HERE, "09_second.py"))
os.environ.setdefault('NMAX', '9000')
s9 = importlib.util.module_from_spec(spec)
import types
src = open(os.path.join(HERE, "09_second.py")).read()
src = src.replace('for lab, r, a, c, d, f, C in d4.ROWS:\n        if lab not in fam:',
                  "SEL = os.environ.get('ROWS', 'R1,R4').split(',')\n    for lab, r, a, c, d, f, C in d4.ROWS:\n        if lab not in fam or lab not in SEL:")
os.environ.setdefault('OUTJSON', 'second_slow.json')
exec(compile(src, '09b', 'exec'), {'__name__': '__main__', '__file__': os.path.join(HERE, '09_second.py')})
