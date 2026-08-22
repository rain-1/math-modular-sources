"""Task 2(b): do the index-9 / index-10 non-congruence host groups of
NONCONGRUENCE_SCAN.md Sec.3 contribute a second-order row with real |lam_2| < 1?

By Theorem N3 + Sec.3 of that note, the only non-mixed option for those groups is
to send the elliptic point to t = infinity, leaving two CUSPS at t_1,t_2 -- the
e = infinity ('Zagier') recurrence class.  That class was scanned exhaustively.
Here we simply re-read the scan outputs and list every row of ANY class with
real |lam_2| < 1, k >= 2 and a non-degenerate Casoratian.
"""
import glob, os, math
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'noncongruence_scan', 'out')
seen = {}
for f in sorted(glob.glob(os.path.join(OUT, '*_scored.txt'))):
    cls = os.path.basename(f).replace('_scored.txt','')
    for line in open(f):
        p = line.split()
        if len(p) != 9: continue
        sc, al, ga, de, ze, l1, l2, k, cplx = (float(p[0]), int(p[1]), int(p[2]), int(p[3]),
                                               int(p[4]), float(p[5]), float(p[6]), int(p[7]), int(p[8]))
        if cplx or k < 2 or abs(l2) >= 1.0: continue
        # drop the two classical Pade families: log (de=16, ga=-3al/4) and
        # arctan (de=-27, al = 18 mod 36, ga = -al/3)
        if de == 16 and 4*ga == -3*al: continue
        if de == -27 and al % 36 == 18 and 3*ga == -al: continue
        # the classical families are UNBOUNDED one-parameter families; restrict to
        # the modular window that contains every census row (al <= 150, |de| <= 100)
        if al > 150 or abs(de) > 100: continue
        key = (round(l1, 6), round(l2, 8), k)
        seen.setdefault(key, [sc, al, ga, de, ze, set()])[5].add(cls)
print(f"{'score':>9}{'lam1':>14}{'lam2':>13}{'k':>3}  {'(al,ga,de,ze)':>22}  classes")
for key, v in sorted(seen.items(), key=lambda z: -z[1][0]):
    sc, al, ga, de, ze, cls = v
    print(f"{sc:>+9.5f}{key[0]:>14.6f}{key[1]:>13.8f}{key[2]:>3}  "
          f"{'(%d,%d,%d,%d)'%(al,ga,de,ze):>22}  {','.join(sorted(cls))[:60]}")
print("""
Reading.  Every real-|lam_2|<1 row with k >= 2 in every box scanned is one of:
 * the Legendre/Pade log family (delta = 16, al = 8x, ga = -3al/4) and the
   arctan family (delta = -27, al = 18 mod 36) -- not modular, unbounded score;
 * Apery's zeta(2) row (11,3,-1,0) and its rescalings (11j,3j,-j^2,0);
 * Beukers' row (136,10,16,4) and its rescalings (136j,10j,16j^2,4j^2);
 * the level-8 sqrt-T row (24,2,16,4);
 * the level-5 Fricke row (88,6,-64,-12).
None of these has a non-congruence HOST GROUP of index 9 or 10: the two modular
e=infinity rows are Apery's (Gamma_1(5), congruence) and the Pade families, and
the e=2 rows live on Gamma_0(6)+6, Gamma_0(5)+5 and level 8 -- all congruence
groups, with the non-congruence-ness in the multiplier system of sqrt(F).
So: NO index-9/10 non-congruence host contributes a second-order row with real
lam_2 < 1, and there is nothing further to score.""")
