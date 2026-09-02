"""19: Landau ceilings and entry margins, one variable vs two, for the three
pullback hosts.  One variable: phi omits the OUTER singularity t2, so
l <= log(16 |t2|)  (Landau; = CDT's 16|s|).  Two variables on a pullback host
u = alpha(x) beta(y) with alpha'(0)=beta'(0)=1 and branch point u_0:
   rho_A rho_B <= |u_0| and rho_A >= e^{l1}/16, rho_B >= e^{l2}/16
   =>  l1 + l2 <= log(256 |u_0|),
plus one-variable Landau bounds from any branch point in x alone or y alone.
"""
from math import log, sqrt

s2 = sqrt(2)
print("%-12s %-26s %8s %6s %9s" % ("row", "one variable", "l<=", "tau", "margin"))
one = [("zeta(3)", "t2=(1+sqrt2)^4=%.3f" % (1+s2)**4, 16*(1+s2)**4, 3),
       ("zeta(2)", "t2=(11+sqrt125)/2=%.3f" % ((11+sqrt(125))/2), 16*(11+sqrt(125))/2, 2),
       ("row E",   "t2=1/4",                                   16*0.25, 2),
       ("s_10",    "t2 = 1/(4(4-sqrt(17))^?) see note",        None, 2)]
for nm, desc, v, tau in one:
    if v is None:
        continue
    print("%-12s %-26s %8.4f %6d %+9.4f" % (nm, desc, log(v), tau, log(v)-tau))
print()
print("%-12s %-30s %8s %8s %10s %10s" %
      ("row", "two variables (pullback)", "l1+l2<=", "min r", "min-margin", "max-margin"))
two = [("zeta(3)", "u0=1/16, alpha=beta=Koebe", 256/16.0, 3),
       ("zeta(2)", "v0=1/4,  alpha=x/(1-x)",    256/4.0,  2),
       ("row E",   "w0=1/4,  beta=y",           256/4.0,  2)]
for nm, desc, v, r in two:
    L = log(v)
    print("%-12s %-30s %8.4f %8d %+10.4f %+10.4f"
          % (nm, desc, L, r, L - r, L - 2*r))
print()
print("extra per-variable Landau bounds (branch points in one variable alone):")
print("   zeta(3): x=1 and y=1 are branch points of 2F1 at u=infinity")
print("            => l_k <= log16 = %.4f  (weaker than the joint bound)" % log(16))
print("   zeta(2): x=1, y=1 likewise => l_k <= log16")
print("   row E  : y=1/4 is a branch point of sqrt(1-4y) => l_2 <= log4 = %.4f;"
      % log(4))
print("            x=1 is a branch point of sqrt(1-4w)   => l_1 <= log16 = %.4f;"
      % log(16))
print("            so max-type entry (needs l_2 > 2) FAILS by 2 - log4 = %+.4f,"
      % (2 - log(4)))
print("            exactly the one-variable Catalan deficit.")
print()
print("min-type entry margins (what a CONVERGENT min-type companion would give):")
for nm, desc, v, r in two:
    print("   %-8s  l1+l2 <= %.4f, need > %d  =>  margin %+0.4f"
          % (nm, log(v), r, log(v)-r))
