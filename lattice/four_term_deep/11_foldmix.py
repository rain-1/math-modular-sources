#!/usr/bin/env python3
"""11_foldmix.py -- fold (connection) constants for MIXED-exponent four-term rows.

Reuses the machinery of 08_fold.py but (i) builds Rc, Sc, Vc from the mixed
normalisation forms of Theorem D3 (see 03_fmix.c / 06_analyse_deep.py), and
(ii) uses the correct local exponent rho_i at each singular point: rho_r at the
rational characteristic root r, rho_p at the two roots of the quadratic factor.

usage:
  python3 11_foldmix.py row  RPN RPD RRN RRD M J1 J2 r a c d f C  [dps1 dps2]
  python3 11_foldmix.py file <analysis.json> <out.txt>   (all CANDIDATE rows)
"""
import sys, os, json, importlib.util
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("f8", os.path.join(HERE, "08_fold.py"))
f8 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f8)
from mpmath import mp, mpf, mpc, nstr, mpmathify

def make_row_mix(cls, row):
    RPN, RPD, RRN, RRD, M, j1, j2 = cls[:7]
    r = cls[7]
    a, c, d, f, C = row
    rp = Fraction(RPN, RPD); rr = Fraction(RRN, RRD)
    g = C * M * M
    assert r != 0 and g % r == 0
    p = g // r; s = a - r
    b = (1 - rr) * r + (1 - rp) * s
    e = -rp * (2 * p + r * s) - rr * r * s
    h = -(1 + 2 * rp + rr) * Fraction(g)
    jj = Fraction(C * j1 * j2)
    assert d == s * r + p
    assert Fraction(j1 + j2, M) == 1 + 2 * rp + rr
    for x in (b, e, h, jj): assert x.denominator == 1
    b, e, h, jj = int(b), int(e), int(h), int(jj)
    Rc = [Fraction(1), Fraction(-a), Fraction(d), Fraction(-g)]
    Sc = [Fraction(1), Fraction(-(a + b)), Fraction(3 * d + e), Fraction(-(5 * g + h))]
    Vc = [Fraction(-c), Fraction(d + e + f), Fraction(-(4 * g + 2 * h + jj))]
    return dict(cls=list(cls), row=list(row), rho=rp, rho_r=rr, rho_p=rp, rroot=r,
                a=a, b=b, c=c, d=d, e=e, f=f, g=g, h=h, j=jj, Rc=Rc, Sc=Sc, Vc=Vc)

def fold_constants_mix(cls, row, dps, verbose=False, frac=mpf('0.4'), dfac=mpf('0.4')):
    old = mp.dps; mp.dps = dps + 25
    try:
        rd = make_row_mix(cls, row)
        rho_r = f8.to_mp(Fraction(rd['rho_r'])); rho_p = f8.to_mp(Fraction(rd['rho_p']))
        rroot = mpf(rd['rroot'])
        P2, P1, P0 = f8.mp_polys(rd)
        ts = f8.sing_points(rd)
        SING = [mpc(0)] + ts
        # the rational root r corresponds to t = 1/r
        tr = mpc(1) / mpc(rroot)
        ir = min(range(3), key=lambda i: abs(ts[i] - tr))
        rhos = [rho_p] * 3; rhos[ir] = rho_r
        Rmin = min(abs(t) for t in ts)
        r0 = mpf('0.45') * Rmin
        NS = int(mp.dps * 2.303 / (-mp.log(mpf('0.45')))) + 200
        ua, ub = f8.taylor_at_zero(rd, NS)
        xis = []; diag = []
        for ic, tc in enumerate(ts):
            rho = rhos[ic]
            dsep = min(abs(tc - s) for s in SING if abs(tc - s) > mp.mpf(10) ** (-(mp.dps // 2)))
            delta = min(dfac * dsep / abs(tc), dfac)
            u = tc / abs(tc); t1 = r0 * u; tstar = tc * (1 - delta); s0 = tstar - tc
            ratio = abs(s0) / dsep
            nloc = min(max(int(mp.dps * 2.303 / (-mp.log(ratio))) + 120, 100), 20000)
            # radial path first; if it runs into another singular point on the ray,
            # retry with a small rotation of the whole path (a detour of the same
            # homotopy class as long as it does not cross a singularity)
            ok = False
            for rot in [mpf(0), mpf('0.35'), mpf('-0.35'), mpf('0.7'), mpf('-0.7')]:
                w = mp.exp(mpc(0, 1) * rot)
                t1r = t1 * w
                pts = ([tstar * w, tstar] if rot != 0 else [tstar])
                try:
                    Av, Ad = f8.series_val(ua, t1r); Bv, Bd = f8.series_val(ub, t1r)
                    _, av, ad, n1 = f8.continue_path(P2, P1, P0, SING, t1r, Av, Ad, 0, pts, frac)
                    _, bv, bd, n2 = f8.continue_path(P2, P1, P0, SING, t1r, Bv, Bd, 1, pts, frac)
                    ok = True; break
                except Exception:
                    continue
            if not ok: raise RuntimeError("continuation failed at t_c")
            kind, w0, w2, part = f8.frob_local(P2, P1, P0, tc, rho, nloc)
            W0, W0p, W1, W1p, PA, PAp = f8.eval_basis(kind, w0, w2, part, s0, rho)
            det = W0 * W1p - W1 * W0p
            c1 = (-av * W0p + ad * W0) / det
            d1 = (-(bv - PA) * W0p + (bd - PAp) * W0) / det
            xis.append(d1 / c1)
            diag.append(dict(tc=tc, kind=kind, rho=rho, c1=c1, d1=d1))
        return [mpmathify(x) for x in ts], [mpmathify(x) for x in xis], diag
    finally:
        mp.dps = old

def fold_checked_mix(cls, row, dps1=110, dps2=170):
    ts1, x1, _ = fold_constants_mix(cls, row, dps1)
    ts2, x2, d2 = fold_constants_mix(cls, row, dps2)
    out = []
    for i, t in enumerate(ts2):
        j = min(range(len(ts1)), key=lambda k: abs(ts1[k] - t))
        out.append((t, x2[i], f8.agree_digits(x1[j], x2[i]), d2[i]))
    return out

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'row':
        v = [int(x) for x in sys.argv[2:15]]
        for i, (t, xi, ag, dg) in enumerate(fold_checked_mix(v[:8], v[8:])):
            print("t_%d = %s %s  rho=%s  agree %d" % (i+1, nstr(mp.re(t),25), nstr(mp.im(t),25),
                                                      nstr(dg['rho'],6), ag))
            print("   Re xi = %s" % nstr(mp.re(xi), 90))
            print("   Im xi = %s" % nstr(mp.im(xi), 90))
    elif cmd == 'file':
        recs = json.load(open(sys.argv[2]))
        fh = open(sys.argv[3], 'w')
        for r in recs:
            if r.get('verdict') != 'CANDIDATE': continue
            lab = ("m" + "_".join(map(str, r['cls'])) + "__" + "_".join(map(str, r['row']))) \
                  if r['kind'] == 'mix' else ("c" + "_".join(map(str, r['cls'])) + "__" + "_".join(map(str, r['row'])))
            try:
                res = fold_checked_mix(r['cls'], r['row']) if r['kind'] == 'mix' \
                      else f8.fold_checked(r['cls'], r['row'], 110, 170)
            except Exception as ex:
                fh.write("# FAIL %s %s\n" % (lab, ex)); fh.flush(); continue
            for i, (t, xi, ag, dg) in enumerate(res):
                fh.write("%s__t%d_re %s\n" % (lab, i+1, nstr(mp.re(xi), max(ag-8, 20))))
                fh.write("%s__t%d_im %s\n" % (lab, i+1, nstr(mp.im(xi), max(ag-8, 20))))
            fh.write("# %s agree %s\n" % (lab, [x[2] for x in res]))
            fh.flush()
        fh.close()
