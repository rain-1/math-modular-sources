"""Host data for CDT_NONCONGREUNCE.md: second-order Apery rows (k = 2), in the
INTEGRAL coordinate x (the one in which a_n in Z).

A row is (al,be,ga,de,ep,ze) with
    (n+1)^2 a_{n+1} = (al n^2 + be n + ga) a_n - (de n^2 + ep n + ze) a_{n-1},
    a_0 = 1, a_1 = ga,
companion b_n same recurrence with b_0 = 0, b_1 = 1 (n >= 1).
Characteristic roots lam_1, lam_2 = roots of L^2 - al L + de, |lam_1| >= |lam_2|.
Singular points of the Picard-Fuchs operator in x:  x_i = 1/lam_i.
score = log(1/|lam_2|) - k    (the elementary Beukers criterion).

'lam' is the geometric denominator of the underlying modular construction
a_n = lam^n [t^n] sqrt(F): the modular coordinate is t = lam*x.  See
geom_denom.py for how CDT's framework absorbs it.
"""
import math
from fractions import Fraction

HOSTS = [
 # name, (al,be,ga,de,ep,ze), k, lam(geometric), field/notes, period, congruence?
 dict(name='Beukers sqrt-Apery (136,10,16,4)', rec=(136,68,10,16,-16,4), k=2, lam=4,
      host='Gamma_0(6)+6, t=(eta_1 eta_6/eta_2 eta_3)^12, F=eta_2^7eta_3^7/(eta_1^5eta_6^5)',
      period='L(Psi,2)=0.1001874492...', cong='NON-congruence (lam=4)'),
 dict(name='level-5 Fricke (88,6,-64,-12)', rec=(88,44,6,-64,64,-12), k=2, lam=2,
      host='Gamma_0(5)+5, t=u/(1+22u+125u^2), u=(eta_5/eta_1)^6, F=E_{2,5}',
      period='0.16430670106434215863... (unidentified)', cong='NON-congruence (lam=2)'),
 dict(name='sqrt-T (24,2,16,4)', rec=(24,12,2,16,-16,4), k=2, lam=1,
      host='level 8, T-row square root', period='L(f_8,2)  [CM weight-3, level 8]',
      cong='congruence (lam=1)'),
 dict(name='sqrt-Domb (20,2,64,16)', rec=(20,10,2,64,-64,16), k=2, lam=1,
      host='level 6, Domb square root', period='L(f_12,2) [CM weight-3, level 12]',
      cong='congruence (lam=1)'),
 dict(name='sqrt-Cooper s_7 (26,2,-27,-6)', rec=(26,13,2,-27,27,-6), k=2, lam=1,
      host='level 7, t=u/(1+13u+49u^2), u=(eta_7/eta_1)^4',
      period='L(f_7,2)  [CM weight-3, level 7]', cong='congruence (lam=1)'),
 dict(name='sqrt-Cooper s_10 (24,2,-256,-60)', rec=(24,12,2,-256,256,-60), k=2, lam=2,
      host='level 10', period='0.31692535921... (unidentified)', cong='NON-congruence (lam=2)'),
 dict(name='sqrt-Cooper s_18 (56,6,768,180)', rec=(56,28,6,768,-768,180), k=2, lam=2,
      host='level 18', period='0.48423755360... (unidentified)', cong='NON-congruence (lam=2)'),
 dict(name='sqrt-AZ(9,3,-27) (72,6,-432,-108)', rec=(72,36,6,-432,432,-108), k=2, lam=4,
      host='AZ(9,3,-27) square root', period='0.14551448201... (unidentified)',
      cong='NON-congruence (lam=4)'),
 # ---- calibration rows (elementary route known / CDT route known) ----------
 dict(name='CAL Apery zeta(2) (11,3,-1,0)', rec=(11,11,3,-1,0,0), k=2, lam=1,
      host='Gamma_1(5) Zagier D', period='zeta(2)/5  [PROVED irrational, elementary]',
      cong='congruence'),
 dict(name='CAL Zagier C = CDT host (10,3,9)', rec=(10,10,3,9,0,0), k=2, lam=1,
      host='Gamma_0(6), CDT L(2,chi_-3)', period='L(2,chi_-3)/2 [PROVED, CDT]',
      cong='congruence'),
 dict(name='CAL Zagier E = Catalan (12,4,32)', rec=(12,12,4,32,0,0), k=2, lam=1,
      host='Gamma_0(8)', period='G = Catalan [OPEN]', cong='congruence'),
]

def roots(rec):
    al, be, ga, de, ep, ze = rec
    D = al*al - 4*de
    if D < 0:
        r = math.sqrt(abs(de)); return (complex(al/2, math.sqrt(-D)/2),
                                        complex(al/2,-math.sqrt(-D)/2), False)
    s = math.sqrt(D)
    r1, r2 = (al+s)/2, (al-s)/2
    if abs(r1) < abs(r2): r1, r2 = r2, r1
    return r1, r2, (D == int(math.isqrt(D))**2)

def enrich(h):
    al, be, ga, de, ep, ze = h['rec']
    l1, l2, sq = roots(h['rec'])
    h['lam1'], h['lam2'] = l1, l2
    h['delta'] = de                                    # = lam_1 lam_2  (norm)
    h['lam2_rational'] = sq                            # discriminant a square?
    h['x2'] = 1.0/l2                                   # outer singular point
    h['x1'] = 1.0/l1
    h['score'] = math.log(1/abs(l2)) - h['k']
    # number-field normalisation |N(lam_2)|^{1/[K:Q]}
    h['lam2_norm'] = abs(l2) if sq else math.sqrt(abs(de))
    return h

for _h in HOSTS: enrich(_h)

if __name__ == '__main__':
    print(f"{'host':<38s} {'lam1':>12s} {'lam2':>12s} {'de':>6s} {'Q?':>3s} "
          f"{'lam2^norm':>9s} {'|x2|':>9s} {'lam':>3s} {'score':>8s}  period")
    for h in HOSTS:
        print(f"{h['name']:<38s} {h['lam1']:>12.6f} {h['lam2']:>12.6f} {h['delta']:>6d} "
              f"{('Q' if h['lam2_rational'] else '-'):>3s} {h['lam2_norm']:>9.4f} "
              f"{abs(h['x2']):>9.5f} {h['lam']:>3d} {h['score']:>+8.4f}  {h['period'][:34]}")
