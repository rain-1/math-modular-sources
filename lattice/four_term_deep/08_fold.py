#!/usr/bin/env python3
"""08_fold.py -- fold (connection) constants of a GENERAL four-term Apery-like row.

A row is (rho; M, j1, j2) + (a, c, d, f, C):

    (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1} + R(n) u_{n-2},  u_0 = 1,
    P = a n^2 + b n + c,   Q = d n^2 + e n + f,   R = g n^2 + h n + j,
    b = (1-rho) a,  e = -2 rho d,  g = C M^2,  h = -C M (j1+j2),  j = C j1 j2.

Its Picard-Fuchs operator (L y = t * [ ... ]) is

    t Rc(t) y'' + Sc(t) y' + Vc(t) y = RHS,
    Rc = 1 - a t + d t^2 - g t^3,
    Sc = 1 - (a+b) t + (3d+e) t^2 - (5g+h) t^3,
    Vc = -c + (d+e+f) t - (4g+2h+j) t^2,

RHS = 0 for A(t) = sum u_n t^n, RHS = 1 for B(t) = sum b_n t^n (b_0=0, b_1=1).

Class identity (holds for every census class): j1 + j2 = (1+3 rho) M, whence

    Sc = Rc + (1-rho) t Rc'                                          (*)

so at a finite singular point t_c (Rc(t_c) = 0, simple) the indicial polynomial
is q21 r (r - rho) with q21 = t_c Rc'(t_c): the local exponents are (0, rho),
exactly as advertised.  (For rho = 0 this degenerates to a double root, giving a
logarithmic / cusp fold; 10_fold.py is that special case.)

FOLD CONSTANT.  Writing near t_c, s = t - t_c,

    rho = 0        : w_0 = 1 + ...,  w_1 = w_0 log s + hh(s)
    rho non-integer: w_0 = 1 + ...,  w_rho = s^rho (1 + ...)

and p(s) a particular solution of the inhomogeneous equation (p(0) = 0), the
connection data are

    A = c0 w_0 + c1 w_*,      B = p + d0 w_0 + d1 w_*,      xi(t_c) = d1 / c1.

xi is independent of the branch of log s / s^rho (a branch change rescales w_*
and hence c1, d1 by the same factor), so it depends only on the homotopy class
of the continuation path.

PATH.  A and B are analytic in |t| < Rmin = min_i |t_i|.  We start at
t1 = r0 * t_c/|t_c| with r0 = 0.45 Rmin (evaluating the t=0 Taylor series
directly there) and continue RADIALLY outwards to t* = (1-delta) t_c.  This is
the canonical "straight-in" path; it never crosses the real axis and never winds
around another singular point.

USAGE
    python3 08_fold.py selftest            # reproduce the K3 row's 3 constants
    python3 08_fold.py census              # all xi=null rows of four_term/out/full.json
    python3 08_fold.py row RN RD M j1 j2 a c d f C
"""
import sys, os, json
from fractions import Fraction
from mpmath import mp, mpf, mpc, log, sqrt, nstr, mpmathify

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, 'out')

# ----------------------------------------------------------------- row -> data
def make_row(cls, row):
    RN, RD, M, j1, j2 = cls
    a, c, d, f, C = row
    rho = Fraction(RN, RD)
    b = (1 - rho) * a
    e = -2 * rho * d
    g = C * M * M
    h = -C * M * (j1 + j2)
    jj = C * j1 * j2
    # class identity check
    assert Fraction(j1 + j2) == (1 + 3 * rho) * M, "class identity j1+j2=(1+3rho)M fails"
    Rc = [Fraction(1), Fraction(-a), Fraction(d), Fraction(-g)]
    Sc = [Fraction(1), -(a + b), 3 * d + e, -(5 * g + h)]
    Vc = [Fraction(-c), d + e + f, -(4 * g + 2 * h + jj)]
    return dict(cls=list(cls), row=list(row), rho=rho,
                a=a, b=b, c=c, d=d, e=e, f=f, g=g, h=h, j=jj,
                Rc=Rc, Sc=Sc, Vc=Vc)

def to_mp(fr):
    return mpf(fr.numerator) / mpf(fr.denominator)

def mp_polys(rd):
    """P2 = t*Rc, P1 = Sc, P0 = Vc as mpc coefficient lists (low -> high)."""
    P2 = [mpc(0)] + [mpc(to_mp(x)) for x in rd['Rc']]
    P1 = [mpc(to_mp(x)) for x in rd['Sc']]
    P0 = [mpc(to_mp(x)) for x in rd['Vc']]
    return P2, P1, P0

def sing_points(rd):
    """the three finite singular points t_i = 1/lam_i, high precision."""
    a, d, g = rd['a'], rd['d'], rd['g']
    coeffs = [mpf(1), mpf(-a), mpf(d), mpf(-g)]          # lam^3 - a lam^2 + d lam - g
    lam = mp.polyroots(coeffs, maxsteps=200, extraprec=4 * mp.prec)
    ts = []
    for L in lam:
        t = mpc(1) / mpc(L)
        # Newton polish on Rc(t) = 1 - a t + d t^2 - g t^3
        for _ in range(60):
            Rv = ((-g * t + d) * t - a) * t + 1
            Rp = (-3 * g * t + 2 * d) * t - a
            dt = Rv / Rp
            t = t - dt
            if abs(dt) < mp.mpf(10) ** (-(mp.dps + 5)) * max(mp.mpf(1), abs(t)):
                break
        ts.append(t)
    return ts

# --------------------------------------------------------------- poly helpers
def shift(poly, t0):
    """Taylor coefficients of poly at t0 (low -> high)."""
    work = [mpc(x) for x in poly]
    out = []
    while work:
        r = work[-1]; newc = [work[-1]]
        for i in range(len(work) - 2, -1, -1):
            r = work[i] + r * t0
            newc.append(r)
        newc.reverse()
        out.append(newc[0])
        work = newc[1:]
    return out

def horner(c, h):
    v = mpc(0)
    for k in range(len(c) - 1, -1, -1):
        v = v * h + c[k]
    return v

def evalser(c, h):
    v = mpc(0)
    for k in range(len(c) - 1, -1, -1):
        v = v * h + c[k]
    d = mpc(0)
    for k in range(len(c) - 1, 0, -1):
        d = d * h + k * c[k]
    return v, d

# -------------------------------------------------- analytic continuation step
def taylor_ordinary(P2, P1, P0, t0, y0, y1, g, nterm):
    A2 = shift(P2, t0); A1 = shift(P1, t0); A0 = shift(P0, t0)
    c = [mpc(y0), mpc(y1)]
    a20 = A2[0]
    for m in range(0, nterm):
        s = mpc(g) if m == 0 else mpc(0)
        for jx in range(1, min(len(A2), m + 3)):
            k = m + 2 - jx
            if 0 <= k < len(c): s -= A2[jx] * c[k] * k * (k - 1)
        for jx in range(0, min(len(A1), m + 2)):
            k = m + 1 - jx
            if 0 <= k < len(c): s -= A1[jx] * c[k] * k
        for jx in range(0, min(len(A0), m + 1)):
            k = m - jx
            if 0 <= k < len(c): s -= A0[jx] * c[k]
        c.append(s / (a20 * (m + 2) * (m + 1)))
    return c

def continue_path(P2, P1, P0, SING, t0, y0, y1, g, pts, frac=mpf('0.4')):
    TOL = mp.mpf(10) ** (-(mp.dps + 15))
    t = mpc(t0); v = mpc(y0); d = mpc(y1)
    nsteps = 0
    for target in pts:
        target = mpc(target)
        while abs(target - t) > TOL:
            R = min(abs(t - s) for s in SING)
            if R == 0:
                raise RuntimeError("continuation hit a singular point")
            step = min(frac * R, abs(target - t))
            h = (target - t) / abs(target - t) * step
            ratio = step / R
            nterm = int(mp.dps * 2.303 / (-mp.log(ratio))) + 60
            nterm = min(max(nterm, 60), 6000)
            c = taylor_ordinary(P2, P1, P0, t, v, d, g, nterm)
            v, d = evalser(c, h)
            t = t + h
            nsteps += 1
            if nsteps > 4000:
                raise RuntimeError("continuation did not terminate")
    return t, v, d, nsteps

# ------------------------------------------------- local Frobenius at t_c
def frob_solve(Q2, Q1, Q0, q21, rho, r, G, c0, nterm):
    """coefficients c_k of s^r sum c_k s^k solving L y = sum G[m] s^{r+m-1}."""
    c = [mpc(c0)]
    for m in range(0, nterm):
        s = mpc(G[m]) if m < len(G) else mpc(0)
        for jx in range(2, min(len(Q2), m + 3)):
            k = m + 2 - jx
            if 0 <= k < len(c): s -= Q2[jx] * c[k] * (r + k) * (r + k - 1)
        for jx in range(1, min(len(Q1), m + 2)):
            k = m + 1 - jx
            if 0 <= k < len(c): s -= Q1[jx] * c[k] * (r + k)
        for jx in range(0, min(len(Q0), m + 1)):
            k = m - jx
            if 0 <= k < len(c): s -= Q0[jx] * c[k]
        c.append(s / (q21 * (r + m + 1) * (r + m + 1 - rho)))
    return c

def frob_local(P2, P1, P0, tc, rho, nterm):
    """returns ('log', w0, hh, part) or ('pow', w0, wr, part)."""
    Q2 = shift(P2, tc); Q1 = shift(P1, tc); Q0 = shift(P0, tc)
    if abs(Q2[0]) > mp.mpf(10) ** (-(mp.dps - 25)) * max(mpf(1), abs(tc)):
        raise RuntimeError("t_c is not a root of t*Rc: %s" % nstr(Q2[0], 8))
    Q2[0] = mpc(0)
    q21 = Q2[1]
    isint = abs(rho - mp.nint(rho)) < mp.mpf(10) ** (-(mp.dps // 2))
    if isint:
        if abs(rho) > mp.mpf('1e-30'):
            raise NotImplementedError("non-zero integer rho not supported")
        w0 = frob_solve(Q2, Q1, Q0, q21, rho, mpf(0), [mpc(0)], 1, nterm)
        R_ = Q2[1:]
        def gq(L, i): return L[i] if 0 <= i < len(L) else mpc(0)
        nd = max(len(Q1), len(Q2))
        D_ = [gq(Q1, jx + 1) - gq(Q2, jx + 2) for jx in range(nd)]
        n = len(w0)
        w0p = [(k + 1) * w0[k + 1] for k in range(n - 1)]
        K = [mpc(0)] * n
        for jx, rj in enumerate(R_):
            for k in range(len(w0p)):
                if jx + k < n: K[jx + k] += 2 * rj * w0p[k]
        for jx, dj in enumerate(D_):
            for k in range(len(w0)):
                if jx + k < n: K[jx + k] += dj * w0[k]
        hh = frob_solve(Q2, Q1, Q0, q21, rho, mpf(0), [-x for x in K], 0, nterm)
        part = frob_solve(Q2, Q1, Q0, q21, rho, mpf(0), [mpc(1)], 0, nterm)
        return ('log', w0, hh, part)
    else:
        w0 = frob_solve(Q2, Q1, Q0, q21, rho, mpf(0), [mpc(0)], 1, nterm)
        wr = frob_solve(Q2, Q1, Q0, q21, rho, rho, [mpc(0)], 1, nterm)
        part = frob_solve(Q2, Q1, Q0, q21, rho, mpf(0), [mpc(1)], 0, nterm)
        return ('pow', w0, wr, part)

def eval_basis(kind, w0, w2, part, s, rho):
    W0, W0p = evalser(w0, s)
    PA, PAp = evalser(part, s)
    if kind == 'log':
        H, Hp = evalser(w2, s)
        L = log(s)
        W1 = W0 * L + H
        W1p = W0p * L + W0 / s + Hp
    else:
        V, Vp = evalser(w2, s)
        sp = mp.exp(rho * log(s))
        W1 = sp * V
        W1p = sp * (Vp + rho * V / s)
    return W0, W0p, W1, W1p, PA, PAp

# ---------------------------------------------------------------- start series
def taylor_at_zero(rd, NS):
    a = mpf(rd['a']); b = to_mp(Fraction(rd['b'])); c = mpf(rd['c'])
    d = mpf(rd['d']); e = to_mp(Fraction(rd['e'])); f = mpf(rd['f'])
    g = mpf(rd['g']); h = mpf(rd['h']); jj = mpf(rd['j'])
    def P(n): return (a * n + b) * n + c
    def Q(n): return (d * n + e) * n + f
    def R(n): return (g * n + h) * n + jj
    aa = [mpf(0), mpf(0), mpf(1)]
    for n in range(0, NS):
        aa.append((P(n) * aa[n + 2] - Q(n) * aa[n + 1] + R(n) * aa[n]) / mpf((n + 1) ** 2))
    bb = [mpf(0), mpf(0), mpf(0), mpf(1)]
    for n in range(1, NS):
        bb.append((P(n) * bb[n + 2] - Q(n) * bb[n + 1] + R(n) * bb[n]) / mpf((n + 1) ** 2))
    return aa[2:], bb[2:]

def series_val(u, t):
    v = mpc(0)
    for k in range(len(u) - 1, -1, -1): v = v * t + u[k]
    dv = mpc(0)
    for k in range(len(u) - 1, 0, -1): dv = dv * t + k * u[k]
    return v, dv

# ---------------------------------------------------------------- main routine
def fold_constants(cls, row, dps, verbose=False, frac=mpf('0.4'), dfac=mpf('0.4')):
    """returns (list of t_i, list of xi_i, diagnostics)"""
    old = mp.dps
    mp.dps = dps + 25
    try:
        rd = make_row(cls, row)
        rho = to_mp(rd['rho'])
        P2, P1, P0 = mp_polys(rd)
        ts = sing_points(rd)
        SING = [mpc(0)] + ts
        Rmin = min(abs(t) for t in ts)
        r0 = mpf('0.45') * Rmin
        NS = int(mp.dps * 2.303 / (-mp.log(mpf('0.45')))) + 200
        ua, ub = taylor_at_zero(rd, NS)
        tail = abs(ua[-1]) * (r0 ** (len(ua) - 1))
        xis = []; diag = []
        for ic, tc in enumerate(ts):
            dsep = min(abs(tc - s) for s in SING if abs(tc - s) > mp.mpf(10) ** (-(mp.dps // 2)))
            delta = min(dfac * dsep / abs(tc), dfac)
            u = tc / abs(tc)
            t1 = r0 * u
            tstar = tc * (1 - delta)
            s0 = tstar - tc
            ratio = abs(s0) / dsep
            nloc = int(mp.dps * 2.303 / (-mp.log(ratio))) + 120
            nloc = min(max(nloc, 100), 20000)
            Av, Ad = series_val(ua, t1)
            Bv, Bd = series_val(ub, t1)
            # residual check of the ODE at t1
            ddA = mpc(0)
            for k in range(len(ua) - 1, 1, -1): ddA = ddA * t1 + k * (k - 1) * ua[k]
            resA = horner(P2, t1) * ddA + horner(P1, t1) * Ad + horner(P0, t1) * Av
            _, av, ad, n1 = continue_path(P2, P1, P0, SING, t1, Av, Ad, 0, [tstar], frac)
            _, bv, bd, n2 = continue_path(P2, P1, P0, SING, t1, Bv, Bd, 1, [tstar], frac)
            kind, w0, w2, part = frob_local(P2, P1, P0, tc, rho, nloc)
            W0, W0p, W1, W1p, PA, PAp = eval_basis(kind, w0, w2, part, s0, rho)
            det = W0 * W1p - W1 * W0p
            c1 = (-av * W0p + ad * W0) / det
            d1 = (-(bv - PA) * W0p + (bd - PAp) * W0) / det
            xi = d1 / c1
            xis.append(xi)
            diag.append(dict(tc=tc, kind=kind, nloc=nloc, nsteps=n1 + n2,
                             resA=abs(resA), tail=tail, c1=c1, d1=d1,
                             locres=abs(w0[-1]) * abs(s0) ** (len(w0) - 1)))
            if verbose:
                print("   t_%d = %s%s  kind=%s nloc=%d steps=%d |c1|=%s"
                      % (ic + 1, nstr(mp.re(tc), 12), nstr(mp.im(tc), 12), kind,
                         nloc, n1 + n2, nstr(abs(c1), 8)))
        res = [mpmathify(x) for x in xis]
        tsc = [mpmathify(x) for x in ts]
        return tsc, res, diag
    finally:
        mp.dps = old

def agree_digits(x, y):
    """number of agreeing decimal digits of two mpc's."""
    with mp.workdps(60):
        dx = abs(x - y); sc = max(abs(x), abs(y))
        if sc == 0: return 10 ** 6
        if dx == 0: return 10 ** 6
        return int(-mp.log10(dx / sc))

def fold_checked(cls, row, dps1=130, dps2=210, verbose=False):
    ts1, x1, d1 = fold_constants(cls, row, dps1, verbose=verbose)
    ts2, x2, d2 = fold_constants(cls, row, dps2, verbose=verbose)
    # match the singular points of the two runs
    out = []
    for i, t in enumerate(ts2):
        j = min(range(len(ts1)), key=lambda k: abs(ts1[k] - t))
        ag = agree_digits(x1[j], x2[i])
        out.append((t, x2[i], ag, d2[i]))
    return out

# --------------------------------------------------------------------- drivers
K3_CLS = [0, 1, 3, 1, 2]
K3_ROW = [11, 4, 37, 3, 3]

def selftest(dps1=130, dps2=210):
    print("=== K3 row self-test: cls %s row %s ===" % (K3_CLS, K3_ROW))
    res = fold_checked(K3_CLS, K3_ROW, dps1, dps2, verbose=True)
    with mp.workdps(dps2):
        Lg = mp.pi / 32 * mp.gamma(mpf(1) / 8) * mp.gamma(mpf(3) / 8) \
             / (mp.gamma(mpf(5) / 8) * mp.gamma(mpf(7) / 8))
        s2 = sqrt(mpf(2))
        targets = {}
        for t, xi, ag, dg in res:
            if abs(mp.im(t)) < mp.mpf('1e-40'):
                pred = (2 * s2 / 3) * Lg
            elif mp.im(t) < 0:
                pred = ((s2 + mpc(0, 1)) / 3) * Lg     # xi((5-i sqrt2)/27) = (sqrt2+i)/3 L ?
            else:
                pred = ((s2 - mpc(0, 1)) / 3) * Lg
            # try both sign conventions
            cands = [pred, mp.conj(pred)]
            best = max(cands, key=lambda p: agree_digits(xi, p))
            targets[(nstr(mp.re(t), 10), nstr(mp.im(t), 10))] = (xi, best,
                                                                 agree_digits(xi, best), ag)
    print()
    ok = True
    for k, (xi, pred, dg, ag) in targets.items():
        print("t = %s %s" % k)
        print("   xi   Re = %s" % nstr(mp.re(xi), 45))
        print("        Im = %s" % nstr(mp.im(xi), 45))
        print("   pred Re = %s" % nstr(mp.re(pred), 45))
        print("        Im = %s" % nstr(mp.im(pred), 45))
        print("   agreement with closed form: %d digits ; dps130-vs-210 self-agreement: %d digits"
              % (dg, ag))
        if dg < 60: ok = False
    print()
    print("SELFTEST %s" % ("PASSED" if ok else "FAILED"))
    return ok

def label(cls, row, i):
    return "c%s__%s__t%d" % ("_".join(str(x) for x in cls),
                             "_".join(str(x) for x in row), i + 1)

def census(dps1=130, dps2=210, outfile=None):
    src = os.path.join(HERE, '..', 'four_term', 'out', 'full.json')
    recs = json.load(open(os.path.abspath(src)))
    todo = [r for r in recs if r.get('xi') is None]
    outfile = outfile or os.path.join(OUT, 'fold_existing.txt')
    f = open(outfile, 'w')
    f.write("# fold constants xi(t_i) of the archimedean-limit-free four-term rows\n")
    f.write("# label  Re(xi)  Im(xi)      (label = c<cls>__<row>__t<i>)\n")
    f.flush()
    nok = nfail = 0
    for r in todo:
        cls, row = r['cls'], r['row']
        print("\n#### cls %s row %s  (rho = %s)" % (cls, row, r['rho']))
        try:
            res = fold_checked(cls, row, dps1, dps2, verbose=True)
        except Exception as ex:
            nfail += 1
            print("   FAILED: %r" % (ex,))
            f.write("# FAILED %s : %r\n" % (label(cls, row, 0), ex)); f.flush()
            continue
        nok += 1
        with mp.workdps(dps2):
            f.write("\n# cls %s row %s rho=%s\n" % (cls, row, r['rho']))
            for i, (t, xi, ag, dg) in enumerate(res):
                nd = min(ag, dps1) - 5
                f.write("# t_%d = %s %s   (agreeing digits: %d)\n"
                        % (i + 1, nstr(mp.re(t), 40), nstr(mp.im(t), 40), ag))
                f.write("%s  %s  %s\n" % (label(cls, row, i),
                                          nstr(mp.re(xi), max(nd, 5)),
                                          nstr(mp.im(xi), max(nd, 5))))
                f.write("%s_re  %s\n" % (label(cls, row, i), nstr(mp.re(xi), max(nd, 5))))
                f.write("%s_im  %s\n" % (label(cls, row, i), nstr(mp.im(xi), max(nd, 5))))
                print("   xi_%d agreeing digits %d : Re=%s Im=%s"
                      % (i + 1, ag, nstr(mp.re(xi), 30), nstr(mp.im(xi), 30)))
            # conjugate pairs -> real sums / differences
            for i in range(len(res)):
                for jx in range(i + 1, len(res)):
                    ti, xii = res[i][0], res[i][1]
                    tj, xij = res[jx][0], res[jx][1]
                    if abs(ti - mp.conj(tj)) < mp.mpf('1e-40') and abs(mp.im(ti)) > mp.mpf('1e-40'):
                        S = xii + xij
                        D = (xii - xij) / mpc(0, 1)
                        ag = min(res[i][2], res[jx][2]); nd = min(ag, dps1) - 5
                        f.write("%s_SUM_t%d_t%d  %s  %s\n"
                                % (label(cls, row, -1).replace('__t0', ''), i + 1, jx + 1,
                                   nstr(mp.re(S), max(nd, 5)), nstr(mp.im(S), max(nd, 5))))
                        f.write("%s_DIFFoverI_t%d_t%d  %s  %s\n"
                                % (label(cls, row, -1).replace('__t0', ''), i + 1, jx + 1,
                                   nstr(mp.re(D), max(nd, 5)), nstr(mp.im(D), max(nd, 5))))
                        print("   xi_%d+xi_%d = %s (+ %s i)" % (i + 1, jx + 1,
                              nstr(mp.re(S), 30), nstr(mp.im(S), 8)))
        f.flush()
    f.close()
    print("\nrows OK: %d   failed: %d   -> %s" % (nok, nfail, outfile))
    return nok, nfail

if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'selftest'
    if cmd == 'selftest':
        d1 = int(sys.argv[2]) if len(sys.argv) > 2 else 130
        d2 = int(sys.argv[3]) if len(sys.argv) > 3 else 210
        ok = selftest(d1, d2)
        sys.exit(0 if ok else 1)
    elif cmd == 'census':
        d1 = int(sys.argv[2]) if len(sys.argv) > 2 else 130
        d2 = int(sys.argv[3]) if len(sys.argv) > 3 else 210
        census(d1, d2)
    elif cmd == 'row':
        v = [int(x) for x in sys.argv[2:12]]
        res = fold_checked(v[:5], v[5:], verbose=True)
        for i, (t, xi, ag, dg) in enumerate(res):
            print("t_%d = %s %s" % (i + 1, nstr(mp.re(t), 30), nstr(mp.im(t), 30)))
            print("   Re xi = %s" % nstr(mp.re(xi), 120))
            print("   Im xi = %s" % nstr(mp.im(xi), 120))
            print("   agreeing digits %d" % ag)
    else:
        print(__doc__)
