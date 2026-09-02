"""10_score.py -- CDT_FINDER scoring (same conventions as lattice/hostscan/11_score.py)
of every imaginary-quadratic candidate produced by this sweep, plus the two reference
points (CDT's own host over Q, and the Gamma_1(5) real-quadratic host) that bracket it.

tau from lattice/cdt_finder/cdt_bound.py with CDT's inventory m=14, u=(1,3[,5]), b_j=2;
ceil = log(256/lam2^norm); entryR = ceil + log(0.62922) - tau;
margin = 14*entryR - (11.845 + log s), s = 1/lam2^norm.
For K imaginary quadratic lam2^norm = |N(lam2)|^{1/2} = |lam2| at the unique (complex)
archimedean place, and the number-field bound is IDENTICAL to the rational one (sec.1).
"""
import math, sys, json
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from cdt_bound import tau_flat, tau_sharp

LOSS = math.log(0.6292232680); BC = 11.845; M = 14

def tau_for(k, m=M):
    u = [max(0, min(m, round((2*j-1)*m/14.0))) for j in range(1, k+1)]
    sm, tf = tau_flat(m, [(uj, 2) for uj in u])
    ec = round(6*m/14.0)
    ts, _ = tau_sharp(m, [1]*ec + [0]*(m-ec))
    return float(tf)+ts, u

def score(l2n, k, m=M):
    s = 1.0/l2n
    T, u = tau_for(k, m)
    ceil_ = math.log(256*s); real = ceil_+LOSS
    return dict(tau=T, ceil=ceil_, entryC=ceil_-T, entryR=real-T,
                margin=m*(real-T)-(BC+math.log(s)), u=u)

r3 = math.sqrt(3)
rows = [
 # name, K, lam2^norm, k, CDT-shape?, live?, period slot, note
 ("CDT / Zagier C   Gamma_0(6) pole 1/2", "Q", 1.0, 2, True, True, "L(2,chi_-3)/2",
  "reference: the theorem"),
 ("Zagier D  Gamma_1(5) pole 2/5",      "Q(sqrt5)", 1.0, 2, True, True, "zeta(2)/5",
  "reference: real quadratic, 2 real places, -1.46 nats tax (GAMMA15_CLOSURE)"),
 ("HYPOTHETICAL imag.quad. host, lam2 a root of unity", "Q(i) or Q(zeta_3)", 1.0, 2, True, True,
  "L(2,psi_4) / L(3,psi_3)", "the prize: NO tax, CDT's own numbers -- does not exist"),
 ("X_0(5) pole at ell_2  (weight 2)",   "Q(i)",      4.0, 3, False, True, "zeta(3)-slot (real)",
  "5-term recurrence: not CDT-shape; lam2 = -4i, N = 16, not a unit"),
 ("X_0(5) pole at ell_2, if k=2 were reached", "Q(i)", 4.0, 2, False, True, "-",
  "hypothetical best case for this geometry"),
 ("X_0(7) pole at ell_3  (weight 1, chi_-7)", "Q(sqrt-3)", 3*r3, 2, False, True, "L(2,chi_-7)-slot (real)",
  "4-term recurrence: not CDT-shape; lam2 = -3 sqrt-3, N = 27"),
 ("X_0(7) pole at ell_3  (weight 2)",   "Q(sqrt-3)", 3*r3, 3, False, True, "zeta(3)-slot (real)",
  "5-term recurrence"),
 ("X_0(9) pole at cusp 1/3 (weight 1)", "Q(sqrt-3)", 3*r3, 2, True, False, "none",
  "CDT-shape but |lam_1| = |lam_2| = 3 sqrt3: no dominant root, no Apery limit"),
 ("X_0(5)/X_0(7)/X_0(9) pole at cusp 0", "Q (coords), lam in imag.quad.", None, 2, None, False, "none",
  "lam_1 = conj(lam_2): equal moduli, dead"),
 ("Zagier E  Gamma_0(8) (Catalan)",     "Q", 4.0, 2, True, True, "G",
  "reference: the nearest Q-rational miss"),
]
out = []
print(f"{'host':52s} {'K':22s} {'l2n':>7s} {'k':>2s} {'tau':>6s} {'ceil':>6s} {'entryC':>7s} {'entryR':>7s} {'margin':>8s}  shape live  period")
for (nm, K, l2n, k, shape, live, per, note) in rows:
    if l2n is None:
        print(f"{nm:52s} {K:22s} {'--':>7s} {k:2d} {'--':>6s} {'--':>6s} {'--':>7s} {'--':>7s} {'--':>8s}  {str(shape):5s} {str(live):5s}  {per}")
        out.append(dict(host=nm, K=K, note=note, dead=True)); continue
    S = score(l2n, k)
    print(f"{nm:52s} {K:22s} {l2n:7.4f} {k:2d} {S['tau']:6.3f} {S['ceil']:6.3f} "
          f"{S['entryC']:+7.3f} {S['entryR']:+7.3f} {S['margin']:+8.3f}  {str(shape):5s} {str(live):5s}  {per}")
    out.append(dict(host=nm, K=K, l2n=l2n, k=k, cdt_shape=shape, live=live, period=per,
                    note=note, **{a: S[a] for a in ('tau','ceil','entryC','entryR','margin')}))
json.dump(out, open('10_scored.json','w'), indent=1)

print()
print("no-tax check (sec.1).  K imaginary quadratic: one archimedean place, d_v = 2 = n_K,")
print("so  m <= (d_v M)/(d_v L - n_K tau) = 2M/(2L-2tau) = M/(L-tau): the rational bound.")
for (lbl, L, num, nK) in [("Q, CDT's own host", 5.081908, 11.845, 1),
                          ("Q(i) host, same geometry (one complex place, d_v=2)", 5.081908, 11.845, 2),
                          ("Q(sqrt5) Gamma_1(5) (two real places, averaged)", 4.90502, 10.81633, 2)]:
    T,_ = tau_for(2)
    den = L - T
    print(f"  {lbl:52s}  L-tau = {den:+.5f}   m <= {num/den:6.3f}   margin(m=14) = {14*den-num:+.4f}")
