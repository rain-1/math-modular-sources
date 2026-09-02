/* 10_cot.gp -- the two structural theorems behind the census.
   INNER (phi = 1):  w(z) = z/(1-z) has Re w = -1/2 for EVERY z on the unit circle
   (and w(1) := zeta(0) = -1/2 by regularisation), Im w(e^{i th}) = (1/2) cot(th/2).
   Hence at every cusp a/c
        Pi = -1/2 L(r,psi) P(r)  +  (i/2) sum_d c_d d^{-r} sum_f psi(f) f^{-r} cot(pi a d f/c),
   the real part is cusp-independent (orientation factor +1, exact Galois
   equivariance across cusps), and the imaginary part is an elementary cotangent sum.
   OUTER (psi = 1):  Pi = sum_f f^{-r} Lam_phi(zeta^f), Lam_phi = N_phi/(1-z^Q), whose
   real part is NOT constant on the unit circle -- whence the factors -3 (A), -11 (D). */
read("lib.gp");
default(realprecision, 50);
FMAX = 400000;
cotsum(s, a, c) = my(r=s[2]); my(psi=s[4]); my(dl=s[6]); my(cl=s[7]); my(tot=0); for(i=1,#dl, my(acc=0); for(f=1,FMAX, my(th=Pi*a*dl[i]*f/c); if(abs(sin(th))>1e-30, acc += pv(psi,f)/f^r * cos(th)/sin(th))); tot += cl[i]*dl[i]^(-r)*acc); tot/2;
tst(nm, a, c) = my(s=srcbyname(nm)); my(b=cperiod(s,a,c)); my(cs=cotsum(s,a,c)); printf("  %-6s cusp %d/%-3d  Im Pi = %s   (i/2) cot-sum = %s   diff = %s\n", nm, a, c, imag(b[1]), cs, abs(imag(b[1])-cs));
print("=== inner rows: imaginary part = the cotangent sum (truncated at f = 400000)");
tst("C",1,3); tst("E",1,4); tst("F",1,6); tst("F",1,3); tst("eps",1,4); tst("alpha",1,3);
print();
print("=== real parts: -1/2 L(r,psi) P(r) at every cusp");
rp(nm) = my(s=srcbyname(nm)); print("  ", nm, " : ", vector(4, j, real(cperiod(s, 1, [1,2,3,4][j])[1])));
rp("C"); rp("E"); rp("F"); rp("alpha");
print();
print("=== outer weight Lam_phi(z) on the unit circle is NOT of constant real part");
WDp = [5,[1,-2,2,-1,0]];
print("  row D, phi = (1,-2,2,-1,0), Q = 5:  Lam(1) = ", wph(WDp,1.0), "   Lam(-1) = ", wph(WDp,-1.0));
print("  factor at cusp 1/2 = [2^-2 Lam(1) + (1-2^-2) Lam(-1)]/Lam(1) = ", (wph(WDp,1.0)/4 + 3*wph(WDp,-1.0)/4)/wph(WDp,1.0));
CH3p = [3,[1,-1,0]];
print("  row A, phi = chi_-3, Q = 3:  Lam(1) = ", wph(CH3p,1.0), "   Lam(-1) = ", wph(CH3p,-1.0));
print("  with the oldform factor (1 - V_2), factor at cusp 1/2 = ", ((wph(CH3p,1.0)/4 + 3*wph(CH3p,-1.0)/4) - wph(CH3p,1.0)/4)/((1-1/4)*wph(CH3p,1.0)));
quit;
