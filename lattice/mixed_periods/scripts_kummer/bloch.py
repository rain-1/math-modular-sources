from mpmath import mp, mpf, mpc, log, exp, pi, polylog, im, re, arg, sqrt, atan, mpmathify, catalan, nstr, pslq, fabs

def BW(z):
    """Bloch-Wigner dilogarithm D(z) = Im Li2(z) + arg(1-z) log|z|"""
    z = mpc(z)
    if abs(z) < mp.eps: return mpf(0)
    if abs(z-1) < mp.eps: return mpf(0)
    return im(polylog(2,z)) + arg(1-z)*log(abs(z))

def Cl2(theta):
    return im(polylog(2, exp(mpc(0,1)*theta)))

def orbit(z):
    """the 6-element anharmonic orbit of z"""
    z=mpc(z)
    return [z, 1-z, 1/z, 1/(1-z), (z-1)/z, z/(z-1)]

def canon(z, pts):
    """find index in pts of an orbit-or-conjugate representative; return (idx, sign) with D(z)=sign*D(pts[idx])"""
    z=mpc(z)
    tol = mpf(10)**(-mp.dps+15)
    for i,p in enumerate(pts):
        for s,w in [(1,z),(-1,mp.conj(z))]:
            o = orbit(w)
            # D(z)=D(1/(1-z))=D((z-1)/z) ;  D(1-z)=D(1/z)=D(z/(z-1))=-D(z)
            sg = [1,-1,-1,1,1,-1]
            for t,q in enumerate(o):
                if abs(q-p) < tol*max(1,abs(p)):
                    return i, s*sg[t]
    return None, None
