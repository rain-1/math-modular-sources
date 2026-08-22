"""Irrationality measures from the CDT arithmetic holonomy bound.

Formula: ICM.tex Remark "withintegrationsremark" (eq. `withintegrations`),
the quantitative refinement of the holonomy bound
(= Theorem `theotrue` with K = Q and S = {v} a singleton):

    m  <=  BC / [ log|phi'(0)| - tau^flat(b) - tau^sharp(e)
                  - (1-gamma)( 2/kappa - (1-gamma)/kappa^2 ) log(1/rho) ]

with
  gamma  = dim_{Q(x)} V(phi,b) / dim_{R(x)} V(phi,b,eta)
         = (number of unconditional "genuine" functions) / m,
  E      = 1 - gamma  (= sum_nu nu*m_nu/m, the "power-system" mean exponent),
  rho    = radius on which the *individual* A_i and B_i converge after pullback,
           i.e. phi(|z|<rho) is inside the disc of convergence of A (the fold).

Writing  Delta0 := m*(log|phi'(0)| - tau) - BC   (the qualitative CDT margin),
the inequality is contradicted as soon as

    Delta0 * kappa^2 - 2*m*L*E*kappa + m*L*E^2 > 0,      L := log(1/rho),

so the resulting measure is the larger root

    kappa_* = E*( m*L + sqrt( m*L*(m*L - Delta0) ) ) / Delta0 ,

and mu_eff(eta) <= kappa_*.  (Delta0 > 0 is exactly the qualitative margin, i.e.
the irrationality proof; kappa_* ~ 2*E*m*L/Delta0 when Delta0 << m*L.)

Calibration: CDT's own L(2,chi_{-3}) measure, ICM Sec. 6.1, kappa = 24781.
"""
import math, json, os
import mpmath as mp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
NCDIR = '/home/ubuntu/code/math-modular-sources/lattice/cdt_noncongruence'

# ----------------------------------------------------------------- kappa
def kappa_star(m, L, E, Delta0):
    """smallest kappa for which the quantitative bound is contradicted."""
    if Delta0 <= 0:
        return None
    mL = m*L
    disc = mL*(mL - Delta0)
    if disc < 0:                      # Delta0 > mL: bound already dead at kappa->E
        return None
    return E*(mL + math.sqrt(disc))/Delta0

# --------------------------------------------------- CDT calibration
def cdt_calibration():
    psi0 = mp.mpf(5448339453535586608)*mp.mpf(10)**9/mp.mpf('8658833407565631122430056127')
    logphi = mp.log(256*psi0)
    tau = mp.mpf(27)/80 + mp.mpf(191)/49
    BC = mp.mpf('11.845')
    m, gamma = 14, mp.mpf(1)/2
    E = 1-gamma
    Delta0 = m*(logphi-tau) - BC
    L = mp.log(11614)
    k = kappa_star(m, float(L), float(E), float(Delta0))
    return dict(logphi=float(logphi), tau=float(tau), entry=float(logphi-tau),
                Delta0=float(Delta0), L=float(L), kappa=k)

# ----------------------------------------------------- conformal data
def lam(w):
    return (mp.jtheta(2, 0, w)/mp.jtheta(3, 0, w))**4

def maxlam(s, N=256):
    """max_{|w|=s} |lambda(w)|  (lambda = 16w -128w^2 + ... , nome w)."""
    s = mp.mpf(s)
    return max(abs(lam(s*mp.e**(2j*mp.pi*mp.mpf(k)/N))) for k in range(N))

def rho_kodaira(r, ratio):
    """phi_r(z) = x_2 lambda(r z).  rho with max_{|z|=rho}|phi_r| = |x_1|,
       ratio = |x_1/x_2| = |lambda_2/lambda_1|."""
    f = lambda s: maxlam(s) - ratio
    s = mp.findroot(f, mp.mpf(ratio)/16)
    return float(s)/r

def rho_koebe(ratio):
    """phi(z) = -4 x_2 z/(1-z)^2 ; max_{|z|=rho}|phi| = 4|x_2| rho/(1-rho)^2 = |x_1|."""
    c = mp.mpf(ratio)/4
    rho = mp.findroot(lambda t: t/(1-t)**2 - c, c)
    return float(rho)

DELTA = json.load(open(os.path.join(NCDIR, 'delta_table.json')))
def delta_of_r(r):
    return DELTA[str(round(r, 2))]

# ------------------------------------------------------------ rows
# (name, lambda_1, lambda_2, k, classical measure, note)
ROWS = [
    ("Beukers (136,10,16,4)", 4*(17+12*math.sqrt(2)), 4*(17-12*math.sqrt(2)), 2, 2,
     50.65365917218, "Beukers 1987 Thm 3; L(Psi,2)=0.10018744..."),
    ("Apery zeta(2) (11,3,-1,0)", (11+5*math.sqrt(5))/2, (11-5*math.sqrt(5))/2, 2, 2,
     11.85078219105, "Apery; Rhin-Viola 5.095"),
    ("Apery zeta(3) (control)", 17+12*math.sqrt(2), 17-12*math.sqrt(2), 3, 3,
     13.41782, "Apery; Rhin-Viola 5.5138"),
]

def tau_flat(m, k):
    """inventory {1} U {theta^j H}: k layers of rate 1, u_j = 1, plus tau^sharp = 0."""
    return k*(1 - 1.0/m/m)

def analyse(name, lam1, lam2, k, R, classical, note):
    # R = holonomic rank; inventory {1} u {theta^j H}_{j<R}, so m = R+1
    x2 = 1.0/abs(lam2); ratio = abs(lam2/lam1)
    out = []
    # (K) univalent
    lp = math.log(4*x2)
    rho = rho_koebe(ratio); L = -math.log(rho)
    for m in (2,):
        tau = tau_flat(m, k)
        D0 = m*(lp-tau) - lp
        E = 1 - 1.0/m
        out.append(("(K) Koebe m=%d" % m, lp, lp, tau, D0, L, E,
                    kappa_star(m, L, E, D0)))
    m = R+1; tau = tau_flat(m, k); E = 1 - 1.0/m
    D0 = m*(lp-tau) - lp
    out.append(("(K) Koebe m=%d" % m, lp, lp, tau, D0, L, E, kappa_star(m, L, E, D0)))
    # (D) Kodaira, m = R+1, optimise over r
    m = R+1; tau = tau_flat(m, k); E = 1 - 1.0/m
    best = None
    for i in range(61):
        r = round(0.20+0.01*i, 2)
        lp = math.log(16*r*x2); BC = lp + delta_of_r(r)
        D0 = m*(lp-tau) - BC
        if D0 <= 0: continue
        rr = rho_kodaira(r, ratio); L = -math.log(rr)
        kk = kappa_star(m, L, E, D0)
        if kk is not None and (best is None or kk < best[-1]):
            best = ("(D) Kodaira m=%d r=%.2f" % (m, r), lp, BC, tau, D0, L, E, kk)
    if best: out.append(best)
    return out

if __name__ == '__main__':
    print(__doc__)
    c = cdt_calibration()
    print("=== calibration: CDT, L(2,chi_-3), ICM Sec.6.1 ===")
    print("  log|phi'(0)| = %.6f   tau = %.6f   entry = %.6f" % (c['logphi'], c['tau'], c['entry']))
    print("  Delta0 (qualitative margin) = %.6f   log(1/rho) = %.6f" % (c['Delta0'], c['L']))
    print("  kappa_*  = %.1f      (CDT print 24781)" % c['kappa'])
    print()
    for row in ROWS:
        name, lam1, lam2, k, R, classical, note = row
        print("=== %s ===" % name)
        print("  lambda_1 = %.6f  lambda_2 = %.8f  k = %d" % (lam1, lam2, k))
        print("  score = log(1/|lam2|) - k = %+.6f" % (math.log(1/abs(lam2))-k))
        print("  classical (Apery argument): mu <= %.5f     [%s]" % (classical, note))
        print("  %-24s %9s %9s %8s %9s %8s %7s %12s" %
              ("architecture", "log|f'(0)|", "BC", "tau", "Delta0", "log1/rho", "E", "kappa_*"))
        for (nm, lp, BC, tau, D0, L, E, kk) in analyse(*row):
            print("  %-24s %9.5f %9.5f %8.5f %9.5f %8.5f %7.4f %12s" %
                  (nm, lp, BC, tau, D0, L, E,
                   ("%.3f" % kk) if kk else "-"))
        print()
