/* 03_showcase.gp -- every number of the showcase page.
   Host Gamma_1(5): x = q prod (1-q^n)^{5(n/5)}, F = sum A_n x^n, A_n = sum_k C(n,k)^2 C(n+k,k).
   Source Phi_new = (1+i phi^5) E_3^{psi_4,1} + c.c., companion B = F D^{-2} Phi_new.
   (1) companion formula vs the inhomogeneous recurrence (exact, n <= 12);
   (2) xi to 60 digits, denominators d_n^2 B_n in Z[phi] (n <= 200), rate phi^{-10};
   (3) K_+ : n A_n phi^{-5n} -> phi^{5/2}/(2 pi 5^{1/4}) (Richardson, 30+ digits);
   (4) the constant term phi_0 of F|_1 S at the cusp 0 (mfslashexpansion) and 5 i phi_0/(2 pi);
   (5) the linear form: (A_n xi - B_n)(-phi^5)^n n = (2/(5 sqrt5))(log n + c_0) + ...;
   (6) the zeta(2)-line companion B_D and its Casoratian constant K_- = 1/(5 sqrt5 K_+).   */
default(parisizemax, 2000000000);
default(realprecision, 1400);
NN = 600;
S5 = Mod(y, y^2 - 5); ph = (1 + S5)/2; ph5 = ph^5; iph5 = 5*ph - 8;   /* phi^{-5} */
if(ph5*iph5 != 1, error("phi^-5 wrong"));
emb(z) = subst(lift(z), y, sqrt(5));
A = vector(NN+1); B = vector(NN+1); BD = vector(NN+1);
A[1] = 1; A[2] = 3; B[1] = 0; B[2] = 2; BD[1] = 0; BD[2] = 1;
{for(n = 1, NN-1,
  A[n+2]  = ((11*n^2+11*n+3)*A[n+1]  + n^2*A[n])/(n+1)^2;
  BD[n+2] = ((11*n^2+11*n+3)*BD[n+1] + n^2*BD[n])/(n+1)^2;
  B[n+2]  = ((11*n^2+11*n+3)*B[n+1]  + n^2*B[n] + 2*(-iph5)^n)/(n+1)^2);}
print("A_n: ", vector(8, i, A[i]));
print("B_n (xi-line, from the recurrence with R_n = 2(-phi^-5)^n): ");
for(i = 1, 9, print("  n=", i-1, "  ", lift(B[i])));
/* (1) companion formula */
NQ = 16;
xq = q*prod(n = 1, NQ, (1 - q^n + O(q^(NQ+1)))^(5*kronecker(n,5)));
Fq = sum(n = 0, NQ, A[n+1]*xq^n);
ps4(n) = if(n%5==0, 0, [1, I, -I, -1][n%5]);
cnew(m) = sumdiv(m, d, (2*real(ps4(m/d)) - 2*ph5*imag(ps4(m/d)))*d^2);
Phinew = sum(m = 1, NQ, cnew(m)*q^m) + O(q^(NQ+1));
print("Phi_new = ", lift(Phinew));
Theta = sum(m = 1, NQ, cnew(m)/m^2*q^m) + O(q^(NQ+1));
peel(G, na) = {my(v = vector(na+1), xp = 1 + O(q^(NQ+1)), c); for(n = 0, na, c = polcoeff(G, n); v[n+1] = c; G -= c*xp; xp *= xq); v;}
Bpeel = peel(Fq*Theta, 12);
print("companion formula - recurrence, n<=12: ", vector(13, i, Bpeel[i] - B[i]));
/* the canonical source: Phi_D = F Dx, companion formula for B_D, and the operator */
Dx = q*deriv(xq, q); PhiD = Fq*Dx; print("Phi_D = F Dx = ", PhiD);
ThD = sum(m = 1, NQ, polcoeff(PhiD, m)/m^2*q^m) + O(q^(NQ+1));
print("B_D formula - recurrence, n<=12: ", vector(13, i, peel(Fq*ThD,12)[i] - BD[i]));
Fpeel = peel(Fq, 12); print("F peeled back - A_n: ", vector(13, i, Fpeel[i] - A[i]));
/* the q-expansion of F and its Eisenstein identification */
print("F(q) = ", Fq + O(q^9));
Fpred = 1 + sum(n = 1, NQ, sumdiv(n, d, 3*real(ps4(d)) + imag(ps4(d)))*q^n);
print("F - [1 + sum_n sum_{d|n}(3 Re psi4(d) + Im psi4(d)) q^n] = ", Fq - Fpred);
/* (2) xi and denominators */
default(realprecision, 420); G5 = znstar(5,1); L2 = lfun(lfuncreate([G5,[1]]), 2); default(realprecision, 1400);
xi400 = ((11+5*sqrt(5))/2)*imag(L2) - real(L2);
rat(n) = emb(B[n+1])/A[n+1];
xi = rat(600);   /* accurate to ~1250 digits; checked against lfun to 400 digits: */
print("xi = ", xi400); print("B_600/A_600 - xi(lfun, 400 digits) = ", xi - xi400);
print("B_n/A_n - xi at n=100,200,400,600: ", [rat(100), rat(200), rat(400), rat(600)] - xi);
print("rate: log|B_n/A_n - xi|/n - log phi^{-10} at n=300,600: ", [log(abs(rat(300)-xi))/300, log(abs(rat(600)-xi))/600] + 10*log((1+sqrt(5))/2));
dn = 1; kmax = 0; bad = 0;
inZphi(z) = my(u = lift(z), a0 = polcoeff(u,0), a1 = polcoeff(u,1)); denominator(2*a0) == 1 && denominator(2*a1) == 1 && (2*a0 - 2*a1) % 2 == 0;
{for(n = 1, 200, dn = lcm(dn, n);
  if(!inZphi(dn^2*B[n+1]), bad++);
  if(inZphi(dn*B[n+1]), kmax++));}
print("d_n^2 B_n in Z[phi] fails for ", bad, " of n<=200;  d_n^1 B_n already integral for ", kmax, " values of n");
/* (3) K_+ by Richardson on n A_n phi^{-5n} */
phr = (1+sqrt(5))/2;
an(n) = n*A[n+1]*phr^(-5*n);
nodes = vector(24, j, NN - 24 + j);
Kp = polinterpolate(vector(24, j, 1/nodes[j]), vector(24, j, an(nodes[j])), 0);
Kpred = phr^(5/2)/(2*Pi*5^(1/4));
print("K_+ (Richardson) = ", Kp); print("phi^{5/2}/(2 pi 5^{1/4}) = ", Kpred, "   diff = ", Kp - Kpred);
/* (4) constant term of F at the cusp 0 */
default(realprecision, 60);
mf = mfinit([5, 1, Mod(2,5)], 4); mfb = mfinit([5,1,Mod(3,5)],4);
E1 = mfbasis(mf); E1b = mfbasis(mfb);
print("dim M_1(5,psi4) = ", #E1, ", dim M_1(5,psibar4) = ", #E1b);
print("basis coefficients: ", mfcoefs(E1[1], 6), "  ", mfcoefs(E1b[1], 6));
/* F = c E1 + c' E1b with rational coefficients: solve from q^1, q^2 */
M = matrix(2,2); M[1,1] = mfcoef(E1[1],1); M[1,2] = mfcoef(E1b[1],1); M[2,1] = mfcoef(E1[1],2); M[2,2] = mfcoef(E1b[1],2);
cc = matsolve(M, [3, 4]~); print("F = ", cc[1], " * E1(psi4) + ", cc[2], " * E1(psibar4)");
tonum(z) = my(u = lift(z)); polcoeff(u,0) + I*polcoeff(u,1);  cn = [tonum(cc[1]), tonum(cc[2])];
chk = vector(10, n, cc[1]*mfcoef(E1[1],n) + cc[2]*mfcoef(E1b[1],n) - polcoeff(Fq, n)); print("check q^1..q^10: ", chk);
S = [0,-1;1,0];
e0 = mfslashexpansion(mf, E1[1], S, 3, 0, &pr); e0b = mfslashexpansion(mfb, E1b[1], S, 3, 0, &prb);
print("E1(psi4)|S: ", e0, "  params ", pr); print("E1(psibar4)|S: ", e0b, "  params ", prb);
phi0 = cn[1]*e0[1] + cn[2]*e0b[1];
print("phi_0 = a_0(F|_1 S) = ", phi0);
print("5 i phi_0/(2 pi) = ", 5*I*phi0/(2*Pi), "   K_+ = ", Kp);
print("phi_0 * 5^{5/4} / (-i phi^{5/2}) = ", phi0*5^(5/4)/(-I*phr^(5/2)));
/* (5) the linear form on the xi-line */
default(realprecision, 1400);
g(n) = (A[n+1]*xi - emb(B[n+1]))*(-phr^5)^n*n;
c1 = 2/(5*sqrt(5));
print("g_n = (A_n xi - B_n)(-phi^5)^n n at n = 50,100,200,400,600: ");
{for(i=1,5, my(n=[50,100,200,400,600][i]); print("   n=",n,"  g_n = ", g(n), "   g_n/log n = ", g(n)/log(n), "  g_n - c1 log n = ", g(n) - c1*log(n)));}
/* fit g_n = c1 log n + c0 + sum_{j=1..K} (p_j log n + q_j)/n^j, unknowns c1', c0, p_j, q_j */
K = 8; nun = 2*K + 2; nds = vector(nun, j, 400 - 6*(nun - j));
Mf = matrix(nun, nun); rhs = vector(nun);
{for(i = 1, nun, my(n = nds[i]); Mf[i,1] = log(n); Mf[i,2] = 1; for(j=1,K, Mf[i,2*j+1] = log(n)/n^j; Mf[i,2*j+2] = 1/n^j); rhs[i] = g(n));}
sol = matsolve(Mf, rhs~);
print("fitted leading coefficient = ", sol[1], "   2/(5 sqrt5) = ", c1, "  diff = ", sol[1] - c1);
print("fitted c_0 = ", sol[2]);
Mf2 = matrix(nun, nun-1); rhs2 = vector(nun);
{for(i = 1, nun, my(n = nds[i]); Mf2[i,1] = 1; for(j=1,K, Mf2[i,2*j] = log(n)/n^j; Mf2[i,2*j+1] = 1/n^j); rhs2[i] = g(n) - c1*log(n));}
sol2 = matsolve(Mf2[1..nun-1,], rhs2[1..nun-1]~);
c0 = sol2[1]; print("c_0 with c1 fixed = ", c0);
default(realprecision, 30);
print("lindep[c0, 1, Euler, log phi, log 5, log 2, log pi, pi, xi, zeta(2)] = ", lindep([c0, 1, Euler, log(phr), log(5), log(2), log(Pi), Pi, xi, zeta(2)], 20));
print("lindep[c0*5sqrt5/2, 1, Euler, log phi, log 5, log 2, log pi] = ", lindep([c0*5*sqrt(5)/2, 1, Euler, log(phr), log(5), log(2), log(Pi)], 20));
/* (6) the zeta(2) line: Casoratian */
default(realprecision, 1400);
wD(n) = (A[n+1]*zeta(2)/5 - BD[n+1])*(-phr^5)^n*n;
nodes2 = vector(24, j, 376 + j);
Km = polinterpolate(vector(24, j, 1/nodes2[j]), vector(24, j, wD(nodes2[j])), 0);
print("K_-(zeta(2) line, B_1 = 1) = ", Km, "   1/(5 sqrt5 K_+) = ", 1/(5*sqrt(5)*Kp), "  2 pi/(5^{5/4} phi^{5/2}) = ", 2*Pi/(5^(5/4)*phr^(5/2)));
print("B_D/A_n - zeta(2)/5 at n=600: ", emb(BD[601])/A[601] - zeta(2)/5);
quit;
