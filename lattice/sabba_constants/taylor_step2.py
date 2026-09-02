"""Task C at high precision: local model with log terms at x=1."""
from mpmath import mp, mpf, nstr, matrix, lu_solve, log as mlog
mp.dps = 80

def taylor_coeffs(x0, ic, a2, N):
    A = x0*(1-x0); B = 1-2*x0; Cc = 2*(1-x0)
    b = [mpf(0)]*(N+5)
    b[0]=ic[0]; b[1]=ic[1]; b[2]=ic[2]/2; b[3]=ic[3]/6
    for m in range(0, N):
        num = ( B*(m+3)*(m+2)*(m+1)*m*b[m+3] + Cc*(m+3)*(m+2)*(m+1)*b[m+3]
                - (m+2)*(m+1)*m*(m-1)*b[m+2] - 2*(m+2)*(m+1)*m*b[m+2] - 4*a2*b[m] )
        b[m+4] = -num/(A*(m+4)*(m+3)*(m+2)*(m+1))
    return b[:N+4]

def evaluate(b, u):
    L=len(b)
    y  = sum(b[n]*u**n for n in range(L-1,-1,-1))
    y1 = sum(n*b[n]*u**(n-1) for n in range(L-1,0,-1))
    y2 = sum(n*(n-1)*b[n]*u**(n-2) for n in range(L-1,1,-1))
    y3 = sum(n*(n-1)*(n-2)*b[n]*u**(n-3) for n in range(L-1,2,-1))
    return [y,y1,y2,y3]

def continue_to(ic, a2, xtarget, N=340):
    x0 = mpf(1)/2; cur = list(ic)
    while True:
        R = min(x0, 1-x0)
        step = min(mpf('0.55')*R, xtarget-x0)
        if step <= 0: break
        cur = evaluate(taylor_coeffs(x0, cur, a2, N), step)
        x0 = x0 + step
        if abs(x0-xtarget) < mpf(10)**(-mp.dps+5): break
    return cur

def limit_at_1(ic, a2, M=9, d0=mpf('0.02')):
    """model y(1-d) = c0 + c1 d + c2 d^2 + sum_{k=3}^{M}(a_k + b_k log d) d^k"""
    basis = lambda d: [mpf(1), d, d**2] + [f for k in range(3,M+1) for f in (d**k, mlog(d)*d**k)]
    ds = [d0/mpf(2)**j for j in range(3+2*(M-2))]
    nb = len(basis(ds[0]))
    ds = ds[:nb]
    A = matrix(nb, nb); rhs = matrix(nb, 1)
    for i,d in enumerate(ds):
        bb = basis(d)
        for j in range(nb): A[i,j] = bb[j]
        rhs[i] = continue_to(ic, a2, 1-d)[0]
    sol = lu_solve(A, rhs)
    return sol[0]

Cref = mpf("0.6287366070989547994355879022526302953834327469083313456011767306080768279905057674386550681519432139")
Dref = mpf("0.5766333897301843923978917497829139257961494352757108398411041918076483544012440202031822641121339437")
for a in (1,-1):
    P = [mpf(1), mpf(0), mpf(4*a), mpf(-8*a)]
    Q = [mpf(0), mpf(1), mpf(0), mpf(-4*a)]
    Pv = limit_at_1(P, 1); Qv = limit_at_1(Q, 1)
    print("a=%+d :  P_a(1^-) = %s"%(a, nstr(Pv,45)))
    print("        Q_a(1^-) = %s"%(nstr(Qv,45)))
    if a==1:  print("        2Q/P - C   = %s"%nstr(2*Qv/Pv - Cref,5))
    else:     print("        P/(2Q) - D = %s"%nstr(Pv/(2*Qv) - Dref,5))
    print("        P/Q  = %s"%nstr(Pv/Qv,40))
