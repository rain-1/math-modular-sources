#!/usr/bin/env python3
"""05_uniformise.py -- conformal size at the cusp 0 of a 4-point orbifold on P^1.

The CDT quantity log|varphi'(0)| is the CUSP SIZE of the host geometry: if Gamma
uniformises the orbifold Omega and the cusp over x=0 is normalised to width 1
(stabiliser tau -> tau+1, q = e^{2 pi i tau}), then x = varphi(q) = c q + O(q^2)
and |varphi'(0)| = |c|.  For P^1 minus {0,s,infty} this is 16|s| (modular lambda);
our post-hypothesis hosts have a FOURTH special point and no closed form, so c is
computed here by numerical Fuchsian uniformisation.

Method
------
* Projective structure  u'' + Q(z) u = 0,
      Q = sum_j [ (1-theta_j^2)/(4 (z-z_j)^2) + B_j/(z-z_j) ],
  theta = 0 at a cusp, 1/2 at an order-2 cone point; z_0 = 0 and z = infinity are
  cusps.  Regularity at infinity (Q ~ 1/(4z^2)) gives two linear conditions on the
  three B_j, leaving ONE real accessory parameter.
* Frobenius basis at 0 (exponents 1/2,1/2, logarithmic):
      u1 = z^{1/2} v1,  u2 = u1 log z + z^{1/2} w,  v1(0)=1, w(0)=0,
  so W := u2/u1 = log z + O(z), tau0 := W/(2 pi i) has tau0 -> tau0+1 around z=0,
  and e^{W} = z (1+O(z)).
* Monodromy by Taylor-step continuation of u'' = -Qu along lassos.
* tau = tau0 + beta uniformises iff every monodromy generator preserves the line
  Im tau0 = h (h = -Im beta).  Two generators give ONE real equation in the
  accessory parameter.  Then q = e^{2 pi i tau} = e^{2 pi i beta} e^{W}, so
  z = e^{-2 pi i beta} q (1+O(q)) and
      |varphi'(0)| = e^{2 pi Im beta} = e^{-2 pi h}.

Validation (`python3 05_uniformise.py test`).
"""
import os, sys, math, cmath, json
import numpy as np

PREC = 30          # digit budget for the Taylor steps (accuracy vs speed)
STEPFRAC = 0.45    # step length as a fraction of the distance to the singular set

# ----------------------------------------------------------------- orbifold
class Orb:
    def __init__(self, pts, acc=0.0):
        """pts = [(z1,theta1),(z2,theta2)] finite non-zero special points."""
        self.pts = [(complex(z), float(th)) for z, th in pts]
        self.sing = [0j] + [z for z, _ in self.pts]
        self.acc = float(acc)
        A = [0.25] + [(1 - th * th) / 4 for _, th in self.pts]
        zs = self.sing
        rhs = 0.25 - sum(A)
        B2 = complex(self.acc)
        B1 = (rhs - B2 * zs[2]) / zs[1]
        B0 = -B1 - B2
        self.A, self.B = A, [B0, B1, B2]

    def rad(self, z0):
        return min(abs(z0 - s) for s in self.sing)

    def rad0(self):
        return min(abs(s) for s in self.sing if s != 0)

    def Qtaylor(self, z0, N):
        k = np.arange(N)
        c = np.zeros(N, dtype=complex)
        for a, b, zj in zip(self.A, self.B, self.sing):
            d = z0 - zj
            p = (-1.0 / d) ** k
            c += (a / (d * d)) * (k + 1) * p + (b / d) * p
        return c

    # ---------------- Frobenius basis at 0 ----------------
    def frob(self, z, NS=None):
        key = (complex(z), NS)
        if getattr(self, '_fc', None) is not None and self._fc[0] == key:
            return self._fc[1]
        R = self.rad0()
        rat = abs(z) / R
        if NS is None:
            NS = int(PREC * 2.303 / max(1e-6, -math.log(rat))) + 30
            NS = min(max(NS, 50), 2000)
        k = np.arange(NS)
        Ra = np.zeros(NS, dtype=complex)
        for a, b, zj in zip(self.A, self.B, self.sing):
            if zj == 0:
                continue
            d = -zj
            p = (-1.0 / d) ** k
            Ra += (a / (d * d)) * (k + 1) * p + (b / d) * p
        B0 = self.B[0]
        av = np.zeros(NS, dtype=complex); bv = np.zeros(NS, dtype=complex)
        av[0] = 1.0
        for m in range(1, NS):
            if m >= 2:
                s2 = np.dot(Ra[:m - 1], av[m - 2::-1])
                t2 = np.dot(Ra[:m - 1], bv[m - 2::-1])
            else:
                s2 = 0j; t2 = 0j
            av[m] = -(B0 * av[m - 1] + s2) / (m * m)
            bv[m] = (-2 * m * av[m] - B0 * bv[m - 1] - t2) / (m * m)
        zp = z ** np.arange(NS)
        v1 = np.dot(av, zp); w = np.dot(bv, zp)
        v1p = np.dot(av[1:] * np.arange(1, NS), zp[:NS - 1])
        wp = np.dot(bv[1:] * np.arange(1, NS), zp[:NS - 1])
        sq = cmath.sqrt(z); lg = cmath.log(z)
        u1 = sq * v1
        u1p = v1 / (2 * sq) + sq * v1p
        u2 = u1 * lg + sq * w
        u2p = u1p * lg + u1 / z + w / (2 * sq) + sq * wp
        res = (complex(u1), complex(u1p), complex(u2), complex(u2p))
        self._fc = (key, res)
        return res

    # ---------------- continuation ----------------
    def step(self, z0, z1, uu):
        R = self.rad(z0); h = z1 - z0
        rat = abs(h) / R
        N = int(PREC * 2.303 / max(1e-6, -math.log(rat))) + 15
        N = min(max(N, 20), 400)
        q = self.Qtaylor(z0, N)
        c = np.zeros(N + 2, dtype=complex)
        c[0] = uu[0]; c[1] = uu[1]
        for n in range(N):
            sv = np.dot(q[:n + 1], c[n::-1])
            c[n + 2] = -sv / ((n + 2) * (n + 1))
        p = h ** np.arange(N + 2)
        u = np.dot(c, p)
        up = np.dot(c[1:] * np.arange(1, N + 2), p[:N + 1])
        return (complex(u), complex(up))

    def follow(self, path, uu, nmax=200000):
        """adaptive walk: from the current point step at most 0.35 * (distance to
        the singular set) towards the next node."""
        n = 0
        for i in range(len(path) - 1):
            z0, z1 = path[i], path[i + 1]
            cur = z0
            while abs(z1 - cur) > 1e-15 * max(1.0, abs(z1)):
                hmax = STEPFRAC * self.rad(cur)
                d = z1 - cur
                h = d if abs(d) <= hmax else d / abs(d) * hmax
                nxt = cur + h
                uu = self.step(cur, nxt, uu)
                cur = nxt
                n += 1
                if n > nmax:
                    raise RuntimeError("path too long")
        return uu

    def lasso(self, idx, zb, nseg=24, frac=0.45):
        """simple lasso from zb to the singular point sing[idx] and back.
        The base point is chosen OFF the real axis so that the straight segment
        to each singularity clears the others -- essential, because conjugating a
        generator by a detour destroys the 'preserves a horizontal line' test."""
        c = self.sing[idx]
        rr = frac * min(abs(c - t) for i, t in enumerate(self.sing) if i != idx)
        d = c - zb
        appr = c - rr * d / abs(d)
        seg = [zb + (appr - zb) * (k / nseg) for k in range(nseg + 1)]
        circ = [c + (appr - c) * cmath.exp(2j * math.pi * k / nseg) for k in range(nseg + 1)]
        return seg + circ[1:] + seg[::-1][1:]

    def clearance(self, zb):
        """min distance from the straight lasso segments to the other singularities."""
        out = []
        for idx in (1, 2):
            c = self.sing[idx]
            for k in range(1, 25):
                z = zb + (c - zb) * (k / 25.0) * 0.9
                out.append(min(abs(z - t) for i, t in enumerate(self.sing) if i != idx))
        return min(out)

    def monodromy(self, idx, zb):
        f = self.frob(zb)
        b11, b12, b21, b22 = f[0], f[2], f[1], f[3]
        det = b11 * b22 - b12 * b21
        out = []
        for (u, up) in [(f[0], f[1]), (f[2], f[3])]:
            v = self.follow(self.lasso(idx, zb), (u, up))
            x = (v[0] * b22 - v[1] * b12) / det
            y = (-v[0] * b21 + v[1] * b11) / det
            out.append((x, y))
        (A, B), (C, D) = out
        return A, B, C, D

    def tau_matrix(self, idx, zb):
        A, B, C, D = self.monodromy(idx, zb)
        tp = 2j * math.pi
        return (D, C / tp, tp * B, A)

    @staticmethod
    def line_height(M):
        p, q, r, s = M
        det = p * s - q * r
        if abs(det) < 1e-200 or not (abs(det) == abs(det)):
            return None
        n = cmath.sqrt(det)
        p, q, r, s = p / n, q / n, r / n, s / n
        tr = p + s
        if abs(r) < 1e-13:
            return None
        if abs(abs(tr) - 2) < 1e-6:
            return ((p - s) / (2 * r)).imag
        return (p / r).imag


def default_base(pts, ang=0.9, frac=0.28):
    r = min(abs(complex(z)) for z, _ in pts)
    return frac * r * cmath.exp(1j * ang)


def heights(pts, acc, zb):
    o = Orb(pts, acc)
    out = []
    for idx in (1, 2):
        M = o.tau_matrix(idx, zb)
        h = Orb.line_height(M)
        det = M[0] * M[3] - M[1] * M[2]
        if abs(det) < 1e-200 or det != det:
            out.append((None, 1e18)); continue
        n = cmath.sqrt(det)
        cc = M[2] / n
        imc = abs(cc.imag) / max(1.0, abs(cc))    # RELATIVE: entries blow up for
        out.append((h, imc))                       # conjugated generators
    return out


def residual(pts, acc, zb):
    (h1, i1), (h2, i2) = heights(pts, acc, zb)
    if h1 is None or h2 is None:
        return None, h1, h2
    if h1 is None or h2 is None:
        return None, h1, h2
    return h1 - h2, h1, h2


def _h1(pts, acc, zb):
    o = Orb(pts, acc)
    M = o.tau_matrix(1, zb)
    return Orb.line_height(M)


def solve_size(pts, zb=None, ngrid=400, nfine=40, span=None, verbose=False,
               prec_scan=14, prec_fine=34):
    """pts = [(z1,th1),(z2,th2)]; cusps at 0 and infinity implicit.
    Returns dict(size, acc, h, imc, resid, ceil16)."""
    global PREC
    if zb is None:
        zb = default_base(pts)
    zb = complex(zb)
    ceil16 = 16.0 * min(abs(complex(z)) for z, _ in pts)
    hlo = -math.log(ceil16) / (2 * math.pi) - 0.35   # size <= 16 min|z| (monotonicity)
    hhi = hlo + 1.9
    if span is None:
        span = 8.0 / (abs(complex(pts[0][0])) * abs(complex(pts[1][0])))
        span = max(span, 1.0)
    old = PREC
    try:
        PREC = prec_scan
        xs = [-span + 2 * span * k / ngrid for k in range(ngrid + 1)]
        h1s = []
        for x in xs:
            try:
                h1s.append(_h1(pts, x, zb))
            except Exception:
                h1s.append(None)
        wins = []
        run = None
        for k, h in enumerate(h1s):
            ok = h is not None and hlo <= h <= hhi
            if ok and run is None:
                run = k
            if not ok and run is not None:
                wins.append((max(run - 4, 0), min(k + 3, ngrid)))
                run = None
        if run is not None:
            wins.append((max(run - 4, 0), ngrid))
        roots = []
        for (i0, i1) in wins:
            a0, a1 = xs[i0], xs[i1]
            ys = [a0 + (a1 - a0) * k / nfine for k in range(nfine + 1)]
            vals = []
            for y in ys:
                try:
                    d, u, v = residual(pts, y, zb)
                except Exception:
                    d = None
                vals.append(d if (d is not None and d == d) else None)
            for k in range(nfine):
                if vals[k] is None or vals[k + 1] is None:
                    continue
                if (vals[k] > 0) == (vals[k + 1] > 0):
                    continue
                lo, hi, flo = ys[k], ys[k + 1], vals[k]
                PREC = prec_fine
                for _ in range(70):
                    mid = 0.5 * (lo + hi)
                    try:
                        w = residual(pts, mid, zb)[0]
                    except Exception:
                        break
                    if w is None or w != w:
                        break
                    if (w > 0) == (flo > 0):
                        lo = mid
                    else:
                        hi = mid
                    if abs(hi - lo) < 1e-15 * max(1.0, abs(lo)):
                        break
                acc = 0.5 * (lo + hi)
                (u1, i1_), (u2, i2_) = heights(pts, acc, zb)
                PREC = prec_scan
                if u1 is None or u2 is None:
                    continue
                roots.append(dict(size=math.exp(-2 * math.pi * u1), acc=acc, h=u1,
                                  imc=max(i1_, i2_), resid=u1 - u2,
                                  logsize=-2 * math.pi * u1, ceil16=ceil16))
    finally:
        PREC = old
    if verbose:
        for r in roots:
            print("   cand acc=%+.10f size=%.8f imc=%.2e" % (r['acc'], r['size'], r['imc']))
    good = [r for r in roots if r['imc'] < 1e-6 and 0 < r['size'] <= ceil16 * 1.000001]
    if not good:
        raise RuntimeError("no Fuchsian root (%d candidates)" % len(roots))
    good.sort(key=lambda r: -r['size'])
    return good[0]


def test():
    print("--- validation of 05_uniformise.py ---")
    print("%-42s %-18s %s" % ("configuration", "computed size", "expected"))
    cases = [([(1.0, 0.0), (-1.0, 0.0)], "cusps {0,1,-1,oo}", "= 4 exactly"),
             ([(1.0, 0.0), (20.0, 0.0)], "cusps {0,1,20,oo}", "< 16"),
             ([(1.0, 0.0), (200.0, 0.0)], "cusps {0,1,200,oo}", "-> 16"),
             ([(1.0, 0.0), (5000.0, 0.0)], "cusps {0,1,5000,oo}", "-> 16"),
             ([(1.0, 0.0), (50.0, 0.5)], "cusp 1, order-2 at 50", "-> 16"),
             ([(1.0, 0.0), (5000.0, 0.5)], "cusp 1, order-2 at 5000", "-> 16"),
             ([(1.0, 0.5), (3000.0, 0.5)], "order-2 at 1 and 3000", "-> 16"),
             ([(2.0, 0.0), (4000.0, 0.0)], "cusps {0,2,4000,oo}", "-> 32"),
             ([(1.0, 0.0), (-1.0, 0.5)], "cusp 1, order-2 at -1", "4 < s < 16")]
    for pts, lab, exp in cases:
        try:
            r = solve_size(pts)
            print("%-42s %-18.10f %s   (acc=%+.8f, |Im c|=%.1e, ceil %.4g)"
                  % (lab, r['size'], exp, r['acc'], r['imc'], r['ceil16']))
        except Exception as e:
            print("%-42s FAILED: %s" % (lab, e))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        test()
