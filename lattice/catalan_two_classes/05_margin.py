"""Adelic margins for the level-16 Catalan host (x = eta(2t)eta(16t)^2/(eta(t)^2 eta(8t)),
F = eta(2t)^4/eta(4t)^2), compared with the level-8 host.

Level-16 host geometry (see CATALAN_TWO_CLASSES.md):
  singular x : -1/4 (fold, lambda_1=4), (-1+-i)/4 (|x|=sqrt2/4), -1/2, infinity
  normaliser descent   sigma(x) = -x/(4x+1) = s x/(x-s) with s = -1/4,
  y = x^2/(x-s) = 4x^2/(4x+1); branch point y=4s=-1; extra point y=-1/2.
  => uniformisation ceiling |phi'(0)| <= 256|s| = 64, log = 4.158883  (same as level 8)
  2-adic slopes in y: pure module 2, doubly-small 1 (measured), conditional 0.
"""
import math, sys
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/adelic_holonomy')
from adelic_bound import adelic, fmt

ceil_ = math.log(256*0.25)                      # 4.158883
real_ = ceil_ + math.log(0.6292232680)          # CDT's realised contour loss, transported
BC    = 11.845 + math.log(0.25)                 # 10.458706

def run(tag, m, cols, e, slopes, lp):
    r = adelic(m, cols, e, ({2:slopes} if slopes else {}), lp, BC)
    fmt(tag, r, m)
    return r

print("BC =", round(BC,6), "  ceiling =", round(ceil_,6), "  realised =", round(real_,6))
print()
print("--- CDT inventory, m=14 (7 pure + 7 conditional): level 8 == level 16 archimedean ---")
cols14=[(1,2),(3,2)]; e14=[0,0,1,0,0,0,0,0,0,1,1,1,1,1]
run("arch only, ceiling            ", 14, cols14, e14, None, ceil_)
run("arch only, realised           ", 14, cols14, e14, None, real_)
run("adelic (pure slope 2), ceiling", 14, cols14, e14, [2]*7+[0]*7, ceil_)
run("adelic (pure slope 2), realised",14, cols14, e14, [2]*7+[0]*7, real_)
print()
print("--- level 16: add the ONE doubly-small function (slope 1 in y, unconditional) ---")
# column shape must stay a step shape: the doubly-small row has the same [1..2n]^2 type
cols15=[(1,2),(3,2)]; e15=e14+[0]
run("m=15, slopes 2^7,1,0^7, ceiling", 15, cols15, e15, [2]*7+[1]+[0]*7, ceil_)
run("m=15, slopes 2^7,1,0^7, realised",15, cols15, e15, [2]*7+[1]+[0]*7, real_)
print()
print("--- how many doubly-small functions would be needed?  (hypothetical d copies) ---")
for d in (1,2,3,4,7,14,21):
    m = 14+d
    e = e14+[0]*d
    run(f"m={m} (d={d} doubly-small, slope 1), ceiling", m, cols14, e, [2]*7+[1]*d+[0]*7, ceil_)
print()
print("--- best-conceivable inventory (u1=u2=m/2, all pure denominator-free) ---")
for m in (14, 20, 30, 50):
    cols=[(m//2,2),(m//2,2)]; e=[0]*m
    run(f"best inventory, m={m}, ceiling", m, cols, e, [2]*(m//2)+[0]*(m//2), ceil_)
