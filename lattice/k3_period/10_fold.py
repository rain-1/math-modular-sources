#!/usr/bin/env python3
"""10_fold.py -- connection constants (complex Apery limits) of the K3 row.

Row:  (n+1)^2 u_{n+1} = (11n^2+11n+4)u_n - (37n^2+3)u_{n-1} + 3(3n-1)(3n-2)u_{n-2}
Operator (after dividing L = theta^2 - tP - ... by t):
    p2 y'' + p1 y' + p0 y = g,   p2 = t(1-t)(1-10t+27t^2),  p1 = p2',
    p0 = -4(1-10t+15t^2);   g = 0 for A = sum u_n t^n,  g = 1 for B = sum b_n t^n
    (companion b_0 = 0, b_1 = 1).
Singular points: t = 0, (5 -+ i sqrt2)/27, 1, infinity.  All finite ones have
exponents (0,0), i.e. a logarithmic (cusp) fold.  The connection constant is
    xi(t_c) = (coefficient of log in B at t_c)/(coefficient of log in A at t_c).
Method: Taylor-series analytic continuation of (y, y') along a polyline from a
point inside |t| < 1/sqrt 27, then matching to the local Frobenius basis at t_c.
"""
import sys
from mpmath import mp, mpf, mpc, log, sqrt, nstr

mp.dps = 230
NLOC = 1400          # terms in the local Frobenius series
TOL  = mp.mpf(10)**(-(mp.dps+15))

# ---------------------------------------------------------------- polynomials
# coefficients low -> high
P2 = [0, 1, -11, 37, -27]
P1 = [1, -22, 111, -108]
P0 = [-4, 40, -60]

def shift(poly, t0):
    """Taylor coefficients of poly at t0 (low -> high), synthetic division."""
    work = [mpc(c) for c in poly]
    out = []
    while work:
        # remainder of work at t0 = value; quotient = work'
        r = work[-1]; newc = [work[-1]]
        for i in range(len(work)-2, -1, -1):
            r = work[i] + r*t0
            newc.append(r)
        newc.reverse()
        out.append(newc[0])
        work = newc[1:]
    return out

SING = [mpc(0), mpc(1), (mpc(5)-mpc(0,1)*sqrt(2))/27, (mpc(5)+mpc(0,1)*sqrt(2))/27]

def dist_sing(t0):
    return min(abs(t0-s) for s in SING)

# --------------------------------------------------------- ordinary-point step
def taylor_ordinary(t0, y0, y1, g, nterm):
    """Taylor coefficients c_k of the solution at the ORDINARY point t0 with
    y(t0)=y0, y'(t0)=y1, RHS = g (constant 0 or 1)."""
    A2 = shift(P2, t0); A1 = shift(P1, t0); A0 = shift(P0, t0)
    c = [mpc(y0), mpc(y1)]
    a20 = A2[0]
    for m in range(0, nterm):
        s = mpc(g) if m == 0 else mpc(0)
        # subtract known terms
        for j in range(1, min(len(A2), m+3)):
            k = m+2-j
            if 0 <= k < len(c): s -= A2[j]*c[k]*k*(k-1)
        for j in range(0, min(len(A1), m+2)):
            k = m+1-j
            if 0 <= k < len(c): s -= A1[j]*c[k]*k
        for j in range(0, min(len(A0), m+1)):
            k = m-j
            if 0 <= k < len(c): s -= A0[j]*c[k]
        c.append(s/(a20*(m+2)*(m+1)))
    return c

def evalser(c, h):
    """value and derivative of sum c_k h^k."""
    v = mpc(0); d = mpc(0); pw = mpc(1)
    for k in range(len(c)-1, 0, -1):
        pass
    # Horner
    v = mpc(0)
    for k in range(len(c)-1, -1, -1):
        v = v*h + c[k]
    d = mpc(0)
    for k in range(len(c)-1, 0, -1):
        d = d*h + k*c[k]
    return v, d

def continue_path(t0, y0, y1, g, pts):
    """analytically continue (y,y') from t0 through the list of points pts."""
    t = mpc(t0); v = mpc(y0); d = mpc(y1)
    for target in pts:
        target = mpc(target)
        while abs(target-t) > TOL:
            R = dist_sing(t)
            step = min(mp.mpf('0.5')*R, abs(target-t))
            h = (target-t)/abs(target-t)*step
            ratio = step/R
            nterm = int(mp.dps*2.303/(-mp.log(ratio))) + 60
            nterm = min(max(nterm, 60), 4000)
            c = taylor_ordinary(t, v, d, g, nterm)
            v, d = evalser(c, h)
            t = t + h
    return t, v, d

# ------------------------------------------------- local Frobenius at t_c (0,0)
def frob_local(tc, nterm=NLOC):
    """returns (w0, hh, part) coefficient lists in s = t-tc.
       w0: analytic solution, w0(0)=1
       hh: w1 = w0*log(s) + hh(s), hh_0 = 0
       part: particular solution of L y = 1, part_0 = 0"""
    Q2 = shift(P2, tc); Q1 = shift(P1, tc); Q0 = shift(P0, tc)
    assert abs(Q2[0]) < mp.mpf(10)**(-(mp.dps-20)), Q2[0]
    Q2[0] = mpc(0)
    q21 = Q2[1]
    # r(s) = q2/s ; d(s) = (q1 - r)/s
    R_ = Q2[1:]
    def gq(L, i): return L[i] if 0 <= i < len(L) else mpc(0)
    nd = max(len(Q1), len(Q2))
    D_ = [gq(Q1, j+1) - gq(Q2, j+2) for j in range(nd)]

    def solve(g, c0):
        c = [mpc(c0)]
        for m in range(0, nterm):
            s = mpc(g[m]) if m < len(g) else mpc(0)
            for j in range(2, min(len(Q2), m+3)):
                k = m+2-j
                if 0 <= k < len(c): s -= Q2[j]*c[k]*k*(k-1)
            for j in range(1, min(len(Q1), m+2)):
                k = m+1-j
                if 0 <= k < len(c): s -= Q1[j]*c[k]*k
            for j in range(0, min(len(Q0), m+1)):
                k = m-j
                if 0 <= k < len(c): s -= Q0[j]*c[k]
            c.append(s/(q21*(m+1)**2))
        return c

    w0 = solve([mpc(0)], 1)
    # K = 2 r w0' + d w0
    n = len(w0)
    w0p = [ (k+1)*w0[k+1] for k in range(n-1) ]
    K = [mpc(0)]*(n)
    for j, rj in enumerate(R_):
        for k in range(len(w0p)):
            if j+k < n: K[j+k] += 2*rj*w0p[k]
    for j, dj in enumerate(D_):
        for k in range(len(w0)):
            if j+k < n: K[j+k] += dj*w0[k]
    hh = solve([-x for x in K], 0)
    part = solve([mpc(1)], 0)
    return w0, hh, part

def logcoeffs(tc, tstar, path, tstart, Aval, Ader, Bval, Bder):
    """continue A,B from tstart to tstar, then match at tc."""
    _, av, ad = continue_path(tstart, Aval, Ader, 0, path)
    _, bv, bd = continue_path(tstart, Bval, Bder, 1, path)
    w0, hh, part = frob_local(tc)
    s = mpc(tstar) - mpc(tc)
    W0, W0p = evalser(w0, s)
    H, Hp = evalser(hh, s)
    PA, PAp = evalser(part, s)
    L = log(s)
    W1  = W0*L + H
    W1p = W0p*L + W0/s + Hp
    det = W0*W1p - W1*W0p
    c0 = ( av*W1p - ad*W1)/det
    c1 = (-av*W0p + ad*W0)/det
    d0 = ( (bv-PA)*W1p - (bd-PAp)*W1)/det
    d1 = (-(bv-PA)*W0p + (bd-PAp)*W0)/det
    return c0, c1, d0, d1

# ------------------------------------------------------------------- start data
def start_series(t1, NS):
    """A, A', B, B' at t1 from the Taylor coefficients at 0."""
    a = [mpf(1)]; b = [mpf(0), mpf(1)]
    # u_{n+1} from (n+1)^2 u_{n+1} = P(n)u_n - Q(n)u_{n-1} + R(n)u_{n-2}
    def P(n): return 11*n*n + 11*n + 4
    def Qf(n): return 37*n*n + 3
    def Rf(n): return 27*n*n - 27*n + 6
    aa = [mpf(0), mpf(0), mpf(1)]      # index shift 2 : u_{-2},u_{-1},u_0
    for n in range(0, NS):
        aa.append((P(n)*aa[n+2] - Qf(n)*aa[n+1] + Rf(n)*aa[n])/mpf((n+1)**2))
    bb = [mpf(0), mpf(0), mpf(0), mpf(1)]   # b_{-2},b_{-1},b_0=0,b_1=1
    for n in range(1, NS):
        bb.append((P(n)*bb[n+2] - Qf(n)*bb[n+1] + Rf(n)*bb[n])/mpf((n+1)**2))
    ua = aa[2:]; ub = bb[2:]
    t = mpc(t1)
    Av = mpc(0); Ad = mpc(0); Bv = mpc(0); Bd = mpc(0)
    for k in range(len(ua)-1, -1, -1):
        Av = Av*t + ua[k]
    for k in range(len(ua)-1, 0, -1):
        Ad = Ad*t + k*ua[k]
    for k in range(len(ub)-1, -1, -1):
        Bv = Bv*t + ub[k]
    for k in range(len(ub)-1, 0, -1):
        Bd = Bd*t + k*ub[k]
    return ua, ub, Av, Ad, Bv, Bd

if __name__ == '__main__':
    NS = 2200
    t1 = mpc('0.06', '-0.14')
    ua, ub, Av, Ad, Bv, Bd = start_series(t1, NS)
    print("u_n :", [int(x) for x in ua[:11]])
    print("b_n :", [nstr(x, 12) for x in ub[:6]])
    print("|t1| =", nstr(abs(t1), 12), " tail ~", nstr(abs(ua[-1]*t1**(len(ua)-1)), 8))
    # sanity: ODE residual at t1
    def polyval(c, t):
        v = mpc(0)
        for k in range(len(c)-1, -1, -1): v = v*t + c[k]
        return v
    # second derivative from series
    def dd(u, t):
        v = mpc(0)
        for k in range(len(u)-1, 1, -1): v = v*t + k*(k-1)*u[k]
        return v
    resA = polyval(P2,t1)*dd(ua,t1) + polyval(P1,t1)*Ad + polyval(P0,t1)*Av
    resB = polyval(P2,t1)*dd(ub,t1) + polyval(P1,t1)*Bd + polyval(P0,t1)*Bv - 1
    print("ODE residual A:", nstr(abs(resA), 6), "  B:", nstr(abs(resB), 6))

    tp = (mpc(5) - mpc(0,1)*sqrt(2))/27      # t_+ = (5 - i sqrt2)/27
    tstar = tp - mpc(0, '0.05')
    c0,c1,d0,d1 = logcoeffs(tp, tstar, [tstar], t1, Av, Ad, Bv, Bd)
    xi = d1/c1
    print()
    print("=== fold at t = (5 - i sqrt2)/27  (I_3) ===")
    print("c1 (log coeff of A) =", nstr(c1, 40))
    print("d1 (log coeff of B) =", nstr(d1, 40))
    print("xi = d1/c1 =")
    print("  Re =", nstr(mp.re(xi), 160))
    print("  Im =", nstr(mp.im(xi), 160))
    with open('out/xi_complex.txt','w') as f:
        f.write(nstr(mp.re(xi), 200)+"\n"+nstr(mp.im(xi), 200)+"\n")
        f.write("c1 "+nstr(mp.re(c1),200)+" "+nstr(mp.im(c1),200)+"\n")
        f.write("d1 "+nstr(mp.re(d1),200)+" "+nstr(mp.im(d1),200)+"\n")

    # --- the real fold at t = 1 (I_6) ------------------------------------
    tc1 = mpc(1)
    path = [mpc('0.15','-0.14'), mpc('0.30','-0.10'), mpc('0.5','-0.05'), mpc('0.75',0), mpc('0.9',0)]
    c0b,c1b,d0b,d1b = logcoeffs(tc1, mpc('0.9'), path, t1, Av, Ad, Bv, Bd)
    xi1 = d1b/c1b
    print()
    print("=== fold at t = 1  (I_6) ===")
    print("c1 =", nstr(c1b, 40))
    print("d1 =", nstr(d1b, 40))
    print("xi_1 = d1/c1 =")
    print("  Re =", nstr(mp.re(xi1), 160))
    print("  Im =", nstr(mp.im(xi1), 160))
    with open('out/xi_one.txt','w') as f:
        f.write(nstr(mp.re(xi1), 200)+"\n"+nstr(mp.im(xi1), 200)+"\n")
        f.write("c1 "+nstr(mp.re(c1b),200)+" "+nstr(mp.im(c1b),200)+"\n")
        f.write("d1 "+nstr(mp.re(d1b),200)+" "+nstr(mp.im(d1b),200)+"\n")
