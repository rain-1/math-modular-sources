/* 01_cuspformula.gp -- independent re-derivation of the cusp-period formula
   (Theorem 3.1 of the paper) and of the numbers  zeta(2)/5, -11 zeta(2)/5 (Zagier D),
   zeta(2)/4, -3 zeta(2)/4 (Zagier A), xi, sigma(xi) (Gamma_1(5), inner line),
   and two imaginary parts (C at 1/3, E at 1/4).  60 digits.
   Conventions:  E_{r+1}^{psi,phi} = c_0 + sum_m ( sum_{e|m} psi(m/e) phi(e) e^r ) q^m,
   source Phi = sum_d c_d E(d tau),  Theta = D^{-r} Phi_{>=1}.
   At q = zeta e^{-t}, zeta = e(a/c):  Theta = rho/t + Pi + O(t).
   Pi  = sum_d c_d d^{-r}   sum_f psi(f) f^{-r}   w_phi(zeta^{df})
   rho = sum_d c_d d^{-r-1} sum_f psi(f) f^{-r-1} n_phi(zeta^{df})
   w_phi(z0) = N(z0)/2 - z0 N'(z0)/Q  and n_phi = N(z0)/Q   if z0^Q = 1,
   w_phi(z0) = N(z0)/(1-z0^Q),        n_phi = 0            otherwise,
   N(z) = sum_{a=1}^{Q} phi(a) z^a.  The f-sum is periodic mod M and is a
   Hurwitz-zeta combination.                                                */
default(realprecision, 70);
e(x) = exp(2*Pi*I*x);
val(ch, j) = ch[2][((j-1) % ch[1]) + 1];           /* ch = [Q, [values]] */
Npol(ch, z) = sum(a=1, ch[1], ch[2][a]*z^a);
Npol1(ch, z) = sum(a=1, ch[1], a*ch[2][a]*z^(a-1));
isroot(ch, z) = abs(z^ch[1] - 1) < 10^(-40);
wph(ch, z) = if(isroot(ch,z), Npol(ch,z)/2 - z*Npol1(ch,z)/ch[1], Npol(ch,z)/(1-z^ch[1]));
nph(ch, z) = if(isroot(ch,z), Npol(ch,z)/ch[1], 0);
/* one block E^{psi,phi}(d tau) at cusp a/c: [Pi, rho] before the d^{-r} factors */
block(psi, phi, r, a, c, d) = {
  my(g = gcd(a*d, c), aa = a*d/g, cc = c/g, M = lcm(cc, psi[1]), SP = 0, SR = 0, z);
  for(j = 1, M, z = e(aa*j/cc);
    SP += val(psi,j) * wph(phi,z) * zetahurwitz(r, j/M);
    SR += val(psi,j) * nph(phi,z) * zetahurwitz(r+1, j/M));
  [SP/M^r, SR/M^(r+1)];
}
period(src, a, c) = { my([nm, r, psi, phi, dl, cl] = src, P = 0, R = 0, b);
  for(i=1, #dl, b = block(psi, phi, r, a, c, dl[i]); P += cl[i]*dl[i]^(-r)*b[1]; R += cl[i]*dl[i]^(-r-1)*b[2]);
  [P, R]; }
TRIV = [1,[1]]; CHM3 = [3,[1,-1,0]]; CHM4 = [4,[1,0,-1,0]];
WD = [5,[1,-2,2,-1,0]];                       /* Re psi_4 - 2 Im psi_4 */
PS4 = [5,[1,I,-I,-1,0]]; PS4B = [5,[1,-I,I,-1,0]];
ph5 = ((1+sqrt(5))/2)^5;
srcD = ["D", 2, TRIV, WD, [1], [1]];
srcA = ["A", 2, TRIV, CHM3, [1,2], [1,-1]];
srcC = ["C", 2, CHM3, TRIV, [1,2], [1,-8]];
srcE = ["E", 2, CHM4, TRIV, [1,2], [1,-8]];
G5 = znstar(5,1); L2 = lfun(lfuncreate([G5,[1]]), 2);        /* L(2,psi_4), psi_4(2)=i */
xi  = ph5*imag(L2) - real(L2);  xis = -real(L2) - imag(L2)/ph5;
report(nm, v, target) = print(nm, " = ", v, "   |diff| = ", abs(v - target));
print("=== Zagier D (outer, phi = (1,-2,2,-1,0)) ===");
p = period(srcD, 0, 1); report("Pi(D,0)   ", p[1], zeta(2)/5); print("   rho = ", p[2]);
p = period(srcD, 1, 2); report("Pi(D,1/2) ", p[1], -11*zeta(2)/5); print("   rho = ", p[2]);
print("   Lambda(1) = ", wph(WD,1), "  Lambda(-1) = ", wph(WD,-1), "  factor = ", (wph(WD,1)/4 + 3/4*wph(WD,-1))/wph(WD,1));
print("=== Zagier A (outer chi_{-3}, oldform 1 - V_2) ===");
p = period(srcA, 0, 1); report("Pi(A,0)   ", p[1], zeta(2)/4);
p = period(srcA, 1, 2); report("Pi(A,1/2) ", p[1], -3*zeta(2)/4);
print("=== Gamma_1(5), inner line Phi_new = (1+i phi^5)E^{psi4,1} + c.c. ===");
srcN  = ["new", 2, [5,[0,0,0,0,0]], TRIV, [1], [1]];   /* placeholder, replaced below */
Pnew(a,c)  = (1+I*ph5)*period(["",2,PS4,TRIV,[1],[1]],a,c) + (1-I*ph5)*period(["",2,PS4B,TRIV,[1],[1]],a,c);
Pnews(a,c) = (1-I/ph5)*period(["",2,PS4,TRIV,[1],[1]],a,c) + (1+I/ph5)*period(["",2,PS4B,TRIV,[1],[1]],a,c);
p = Pnew(0,1);  report("Pi(new,0)   ", p[1], xi);  print("   rho(new,0)   = ", p[2]);
p = Pnew(1,2);  report("Pi(new,1/2) ", p[1], xi);  print("   rho(new,1/2) = ", p[2]);
p = Pnews(0,1); report("Pi(new',0)  ", p[1], xis); print("   rho(new',0)  = ", p[2]);
p = Pnews(1,2); report("Pi(new',1/2)", p[1], xis); print("   rho(new',1/2)= ", p[2]);
print("xi = ", xi); print("sigma(xi) = ", xis);
print("=== imaginary parts: C at 1/3, E at 1/4 ===");
G3 = znstar(3,1); L2m3 = lfun(lfuncreate([G3,[1]]),2);
p = period(srcC, 1, 3); report("Pi(C,1/3) ", p[1], L2m3/2 + 2*Pi^2/(9*sqrt(3))*I);
p = period(srcE, 1, 4); report("Pi(E,1/4) ", p[1], Catalan/2 + Pi^2/16*I);
print("=== direct Abel cross-check of -11 zeta(2)/5 (sum to m = 40000, one Richardson step) ===");
cD(m) = sumdiv(m, d, val(WD, d)*d^2);
S(t) = sum(m=1, 40000, cD(m)/m^2*(-1)^m*exp(-m*t));
ab = 2*S(0.00125) - S(0.0025); print("Abel-Richardson = ", ab, "   |diff| = ", abs(ab + 11*zeta(2)/5));
quit;
