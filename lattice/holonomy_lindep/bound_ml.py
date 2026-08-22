"""Task A: the CDT arithmetic holonomy bound applied to MULTI-PERIOD rows.

SET-UP.  A row of holonomic rank R (= theta-order of its operator L, MUM at 0)
with r-1 companions B^(1..r-1) (Lemma 1.1 of MULTISLOPE_PROGRAM.md), all of
denominator type [1,..,n]^k, and archimedean limits xi_1,..,xi_{r-1}.
Hypothesis to be contradicted:  a_0 + sum_i a_i xi_i = 0,  a_i in Q (or in K).
It produces exactly ONE conditional generator

    H = sum_i a_i B^(i) + a_0 A ,     h_n = sum_i a_i (b^(i)_n - xi_i a_n),

*independently of how many periods the relation involves*.  Its theta-orbit is
{H, theta H, ..., theta^{R-1} H} (theta^R H is Q(x)-dependent since L(H) is a
polynomial), so the inventory is {1} u {theta^j H} and

    m = R + 1 ,   tau = tau^flat = k(1 - 1/m^2) ,   tau^sharp = 0 .

FOLD-REGULARITY (Prop. 2.1 of HOLONOMY_LINDEP.md).  |h_n| = O(|lambda_2|^n)
*only if the dominant characteristic root lambda_1 is SIMPLE*.  If lambda_1 has
multiplicity mu, the local solution space at x_1 = 1/lambda_1 carries a Jordan
block of size mu and killing the whole log-tower is a codimension-mu condition;
one rational relation supplies ONE condition.  Verified: AESZ 207 (mu=1) folds
(|h_n| = c 53248^n/n^4); Sym^3(Zagier E) (mu=3) and X_1(5) Sym^2 (mu=2) do NOT.

GEOMETRY.  Omega = C minus the singularities that remain after the fold is
deleted.  If >= 2 finite points remain, Omega is hyperbolic and
|phi'(0)| <= r_Omega(0) =: C for EVERY admissible phi (Schwarz-Pick): a genuine
ceiling.  If exactly 1 remains, Omega is not hyperbolic and the ceiling is
Kodaira's 16|x_2| along the modular-lambda family phi_r = x_2 lambda(rz).

READOUT.  entry = log|phi'(0)| - tau ;  margin = m*entry - BC ;
deficit = margin/(m-1)  [nats per function, the CDT_NONCONGRUENCE.md unit].
HEADROOM = log C - (3/4)k  is a *necessary* condition for the architecture ever
to work with any inventory: the margin (m-1)log C - m tau is positive for large m
only if log C > tau_infinity, and tau_infinity >= (3/4)k for the best
inventory that CDT_FINDER.md Sec.4 scenario B allows (u_j = m/2).
"""
import json, math, os, sys
import mpmath as mp
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from conformal import r_two_punctures, r_upper_multi, r_lower_onecut
mp.mp.dps = 30

NCDIR = '/home/ubuntu/code/math-modular-sources/lattice/cdt_noncongruence'
DELTA = json.load(open(os.path.join(NCDIR, 'delta_table.json')))
DELTA.update(json.load(open('delta_ext.json')))
RS = sorted(float(k) for k, v in DELTA.items() if v > -1e30)

def tau_flat(m, k, u=1):
    return k*(1 - float(u*u)/(m*m))

S2, S5, S17 = mp.sqrt(2), mp.sqrt(5), mp.sqrt(17)
z207_1 = mp.mpf(1)/53248
z207_p = (mp.mpf(349)+85*S17)/mp.mpf(2)**17
z207_m = (mp.mpf(349)-85*S17)/mp.mpf(2)**17
c243 = mp.polyroots([1, 57, -289, -1])          # -61.6848, -0.0034578(fold), 4.68830
c313 = mp.polyroots([1, mp.mpf(37)/6, -mp.mpf(1)/48])
i104 = mp.polyroots([1, 0, mp.mpf(1)/72])

# name, R, ncomp, k, lam1, lam2, simple dominant?, prk, fold, remaining, note
ROWS = [
 ("CALIB Zagier C (10,3,9) = CDT host, Gamma_0(6)", 2,1,2, 9.0,1.0, True, 2,
   mp.mpf(1)/9, [mp.mpf(1)], "1,zeta(2),L(2,chi_-3): 3 cusps -> 3 Eisenstein classes; CDT's theorem"),
 ("CALIB Apery zeta(2) (11,3,-1,0)", 2,1,2, float((11+5*mp.sqrt(5))/2), float((5*mp.sqrt(5)-11)/2), True, 1,
   2/(11+5*S5), [-2/(5*S5-11)], "single period (Gamma_1(5) carries no 2nd real class)"),
 ("CALIB Beukers (136,10,16,4)", 2,1,2, float(4*(17+12*mp.sqrt(2))), float(4*(17-12*mp.sqrt(2))), True, 1,
   1/(4*(17+12*S2)), [1/(4*(17-12*S2))], "single period; Gamma_0(6)+6 has 2 cusp orbits, c<=1"),
 ("CALIB Apery zeta(3)", 3,1,3, float(17+12*mp.sqrt(2)), float(17-12*mp.sqrt(2)), True, 1,
   1/(17+12*S2), [1/(17-12*S2)], "single period"),
 ("AESZ 207 (4.4.38) rank4, 3 companions", 4,3,4, 89531.389206720147, 53248.0, True, 3,
   z207_m, [z207_1, z207_p], "prk=3 at BOTH places; no Q- or Q(sqrt17)-relation below 10^337"),
 ("AESZ 243 (4.5.46) rank4, 4 companions", 4,4,4, 289.197085465, 3.0, True, 4,
   c243[1], [mp.mpf(-1)/3, c243[0], c243[2]], "prk=4; best literal score in the MUM sweep"),
 ("AESZ 388 (4.3.18) rank4, 2 companions", 4,2,4, 1156.0, 4.0, True, 2,
   mp.mpf(1)/1156, [mp.mpf(1)/4], "prk=2 to height 10^82; companion k = (3,4), H has k=4"),
 ("AESZ 313 (4.5.82) rank4, 4 companions", 4,4,4, 296.162, 6.0, True, 4,
   c313[1], [mp.mpf(-1)/6, mp.mpf(1), c313[0]], "prk=4"),
 ("AESZ 34 (4.3.1) rank4, 2 companions", 4,2,4, 25.0, 9.0, True, 2,
   mp.mpf(1)/25, [mp.mpf(1)/9, mp.mpf(1)], "prk=2"),
 ("AESZ 104 (4.8.2) rank4, 7 companions", 4,7,4, 72.0, 9.0, True, 7,
   mp.mpf(1)/72, [mp.mpf(-1), mp.mpf(1)/8, mp.mpf(-1)/9, i104[0], i104[1]], "prk=7 (the largest in the corpus)"),
 ("AESZ 380 (4.5.111) rank4, 4 companions", 4,4,4, 108.0, 10.0, True, 4,
   mp.mpf(1)/108, [mp.mpf(1)/10, mp.mpf(-1)/4], "prk=4"),
 ("AESZ 410 (4.3.28) rank4, 2 companions", 4,2,4, 81.0, 32.0, True, 2,
   mp.mpf(1)/81, [mp.mpf(1)/32], "prk=2"),
 ("AESZ 392 (4.3.22) rank4, 2 companions, k=3", 4,2,3, 676.0, 108.0, True, 2,
   mp.mpf(1)/676, [mp.mpf(-1)/108], "prk=2; k=3 for both companions (rare)"),
 ("zeta(5) level 16 (joint operator rank 11)", 11,1,5, float(2+4*mp.sqrt(2)), 4.0, True, 1,
   1/(2+4*S2), [mp.mpf(-1)/2, mp.mpf(-1)/4, 1/(2-4*S2)], "higher companions have NO k<=14: inadmissible"),
 ("zeta(7) level-12 parent (level-60 programme)", None,3,7, 10.0, 6.0, True, 1,
   mp.mpf(1)/10, [mp.mpf(-1)/6], "3 sources but all periods in Q.zeta(7): prk=1"),
 ("beta(4)@24 = Sym^3(Zagier E) rank4, 5 comp.", 4,5,4, 8.0, 4.0, False, 3,
   mp.mpf(1)/8, [mp.mpf(1)/4], "lambda_1 TRIPLE: no fold. prk=3 (2 genuine Q-relations)"),
 ("X_1(5) Sym^2 rank3, 3 companions", 3,3,3, float((11+5*mp.sqrt(5))/2), float((5*mp.sqrt(5)-11)/2), False, None,
   2/(11+5*S5), [-2/(5*S5-11)], "lambda_1 DOUBLE: no fold (measured here)"),
 ("zeta(7) level 24 = Sym^6 rank7, 2 comp.", 7,2,7, float(4+2*mp.sqrt(2)), float(2+2*mp.sqrt(2)), False, None,
   1/(4+2*S2), [1/(-2-2*S2), mp.mpf(1), 1/(4-2*S2), 1/(-2+2*S2)], "roots of multiplicity 7: no fold"),
]

def geom(rest):
    lo, how = r_lower_onecut(rest, 0)
    if len(rest) == 1:
        return float(lo), how, float(16*abs(rest[0])), 'Kodaira 16|x2|', True
    return float(lo), how, float(r_upper_multi(rest, 0)), 'hyperbolic', False

def run():
    out = []
    for (name,R,nc,k,l1,l2,simple,prk,fold,rest,note) in ROWS:
        rho, how, C, ckind, kod = geom(rest)
        score = math.log(1.0/abs(l2)) - k
        head = math.log(C) - 0.75*k
        m = (R+1) if R else None
        rec = dict(name=name,R=R,nc=nc,k=k,l1=l1,l2=l2,simple=simple,prk=prk,
                   score=score,rho=rho,how=how,C=C,ckind=ckind,kod=kod,head=head,
                   m=m,note=note,fold=float(abs(fold)))
        if m:
            tau = tau_flat(m,k); rec['tau']=tau
            rec['entryK'] = math.log(rho)-tau
            rec['marK'] = m*rec['entryK'] - math.log(rho)
            rec['defK'] = rec['marK']/(m-1)
            rec['entryC'] = math.log(C)-tau
            rec['marC'] = m*rec['entryC'] - math.log(C)     # upper bound
            rec['defC'] = rec['marC']/(m-1)
            if kod:
                best=None
                for r in RS:
                    lp = math.log(16*r/abs(l2)); BC = lp + DELTA[str(r)]
                    d = (m*(lp-tau)-BC)/(m-1)
                    if best is None or d>best[0]: best=(d,r,lp,BC)
                rec['defD'], rec['rD'] = best[0], best[1]
        out.append(rec)
    return out

if __name__ == '__main__':
    print(__doc__)
    res = run()
    print("%-46s %3s %3s %2s %4s %10s %8s %8s %8s %8s %8s" %
          ("row","R","nc","k","prk","score","log rho","log C","defK","defC","HEADROOM"))
    for r in res:
        print("%-46s %3s %3s %2d %4s %+10.4f %8.4f %8.4f %8s %8s %+8.4f" %
              (r['name'][:46], r['R'], r['nc'], r['k'], r['prk'], r['score'],
               math.log(r['rho']), math.log(r['C']),
               ("%.4f"%r['defK']) if r.get('defK') is not None else "-",
               ("%.4f"%r['defC']) if r.get('defC') is not None else "-",
               r['head']))
    print()
    print("fold-regular?  ", {r['name'][:24]: r['simple'] for r in res if not r['simple']})
    print()
    print("Kodaira (D) deficits where exactly one puncture remains:")
    for r in res:
        if r.get('defD') is not None:
            print("   %-46s r*=%.2f  deficit %+8.4f" % (r['name'][:46], r['rD'], r['defD']))
    print()
    print("RANKING by best available deficit (fold-regular rows only):")
    fr = [r for r in res if r['simple'] and r.get('defK') is not None]
    fr.sort(key=lambda r: -max(r.get('defD', -99), r['defC']))
    for r in fr:
        print("   %+8.4f   %-52s  [prk %s]" % (max(r.get('defD',-99), r['defC']), r['name'][:52], r['prk']))
