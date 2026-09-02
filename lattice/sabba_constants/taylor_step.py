"""Analytic continuation of  L_a[y] = x(1-x) y'''' + 2(1-x) y''' - 4 a^2 y = 0
by Taylor stepping, from x=1/2 towards x=1.  Task C: verify GPT-5.6's claim
   C = 2 Q_1(1^-)/P_1(1^-),  D = P_{-1}(1^-)/(2 Q_{-1}(1^-))
with P_a: (y,y',y'',y''')(1/2) = (1,0,4a,-8a),  Q_a: (0,1,0,-4a).
"""
from mpmath import mp, mpf, nstr, matrix
mp.dps = 60

def taylor_coeffs(x0, ic, a2, N):
    """ic = (y,y',y'',y''') at x0. returns b[0..N] with y = sum b_n (x-x0)^n."""
    A = x0*(1-x0); B = 1-2*x0; Cc = 2*(1-x0)
    b = [mpf(0)]*(N+5)
    b[0]=ic[0]; b[1]=ic[1]; b[2]=ic[2]/2; b[3]=ic[3]/6
    for m in range(0, N):
        num = ( B*(m+3)*(m+2)*(m+1)*m*b[m+3] + Cc*(m+3)*(m+2)*(m+1)*b[m+3]
                - (m+2)*(m+1)*m*(m-1)*b[m+2] - 2*(m+2)*(m+1)*m*b[m+2] - 4*a2*b[m] )
        b[m+4] = -num/(A*(m+4)*(m+3)*(m+2)*(m+1))
    return b[:N+4]

def evaluate(b, u):
    """value and 3 derivatives at x0+u"""
    y = sum(b[n]*u**n for n in range(len(b)-1,-1,-1))
    y1 = sum(n*b[n]*u**(n-1) for n in range(len(b)-1,0,-1))
    y2 = sum(n*(n-1)*b[n]*u**(n-2) for n in range(len(b)-1,1,-1))
    y3 = sum(n*(n-1)*(n-2)*b[n]*u**(n-3) for n in range(len(b)-1,2,-1))
    return [y,y1,y2,y3]

def continue_to(ic, a2, xtarget, N=260):
    x0 = mpf(1)/2; cur = list(ic)
    while True:
        # radius of convergence at x0 is min(x0, 1-x0); step 60% of it toward 1
        R = min(x0, 1-x0)
        step = min(mpf('0.6')*R, xtarget-x0)
        if step <= 0: break
        b = taylor_coeffs(x0, cur, a2, N)
        cur = evaluate(b, step)
        x0 = x0 + step
        if abs(x0-xtarget) < mpf(10)**(-mp.dps+5): break
    return cur

def limit_at_1(ic, a2, deltas=(mpf('1e-4'),mpf('5e-5'),mpf('25e-6'),mpf('125e-7'))):
    vals=[]
    for d in deltas:
        v = continue_to(ic, a2, 1-d)[0]
        vals.append((d,v))
    # Richardson: v(d) = c0 + c1 d + c2 d^2 + O(d^3 log d);  Neville in d -> 0
    K=len(vals); h=[v[0] for v in vals]; T=[[mpf(0)]*K for _ in range(K)]
    for i in range(K): T[i][0]=vals[i][1]
    for j in range(1,K):
        for i in range(K-j):
            T[i][j] = (T[i+1][j-1]*(0-h[i]) - T[i][j-1]*(0-h[i+j]))/(h[i+j]-h[i])
    return T[0][K-1], [nstr(v[1],20) for v in vals]

Cref = mpf("0.628736607098954799435587902252630295383432746908331345601176730608076827990506")
Dref = mpf("0.576633389730184392397891749782913925796149435275710839841104191807648354401244")
for a in (1,-1):
    a2 = 1
    P = [mpf(1), mpf(0), mpf(4*a), mpf(-8*a)]
    Q = [mpf(0), mpf(1), mpf(0), mpf(-4*a)]
    Pv,_ = limit_at_1(P, a2); Qv,rows = limit_at_1(Q, a2)
    print("a=%+d :  P_a(1^-) = %s"%(a, nstr(Pv,30)))
    print("        Q_a(1^-) = %s"%(nstr(Qv,30)))
    print("        2Q/P     = %s"%(nstr(2*Qv/Pv,30)))
    print("        P/(2Q)   = %s"%(nstr(Pv/(2*Qv),30)))
    if a==1:  print("        2Q/P - C = %s"%nstr(2*Qv/Pv - Cref,5))
    else:     print("        P/(2Q) - D = %s"%nstr(Pv/(2*Qv) - Dref,5))
