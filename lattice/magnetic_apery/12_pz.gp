default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 80);
NQ = 405; LIM = 400;
E4 = 1 + 240*sum(n=1, NQ-1, sigma(n,3)*q^n) + O(q^NQ);
E6 = 1 - 504*sum(n=1, NQ-1, sigma(n,5)*q^n) + O(q^NQ);
DE = (E4^3 - E6^2)/1728;
print("### Task 5: Pasol-Zudilin's level-1 magnetic weight-4 forms");
Pa = DE/E4^2;
Pb = E4*DE/E6^2;
print("F_4a = Delta/E_4^2 : c(1..10) = ", vector(10,j,polcoeff(Pa,j)));
print("F_4b = E_4 Delta/E_6^2 : c(1..10) = ", vector(10,j,polcoeff(Pb,j)));
Xa = Dinv(Pa, LIM); Xb = Dinv(Pb, LIM);
print("Xi_a integral? first failure m = ", intfail(Xa, LIM), "   (0 = integral to m=400)");
print("Xi_b integral? first failure m = ", intfail(Xb, LIM));
print("magnetic (m | c(m)) : F_4a first failure m = ", magfail(Pa, LIM), " ; F_4b first failure m = ", magfail(Pb, LIM));
print("");
print("growth: |c(m)|^(1/m) at m=400 : a: ", abs(polcoeff(Pa,400)*1.0)^(1/400), "   b: ", abs(polcoeff(Pb,400)*1.0)^(1/400));
print("  (exp(pi*sqrt3) = ", exp(Pi*sqrt(3)), " for the disc -3 pole at rho ; exp(2pi) = ", exp(2*Pi), " for the disc -4 pole at i)");
print("");
\\ Lambda(Phi,s) = int_0^inf y^{s-1} Phi(iy) dy ; Phi|_4 S = Phi (level 1) so Lambda(s)=Lambda(4-s)
\\ split at y=1 : Lambda(s) = sum_m c(m) (2 pi m)^{-s} Gamma(s,2 pi m) + same with s -> 4-s
Lam(cv, s, mm) = sum(m=1, mm, cv[m]*(2*Pi*m)^(-s)*incgam(s, 2*Pi*m)) + sum(m=1, mm, cv[m]*(2*Pi*m)^(s-4)*incgam(4-s, 2*Pi*m));
cva = vector(LIM, m, polcoeff(Pa,m));
{
my(L3, L2, L1, xi);
L3 = Lam(cva, 3, LIM); L2 = Lam(cva, 2, LIM); L1 = Lam(cva, 1, LIM);
print("F_4a : Lambda(2) = ", L2);
print("F_4a : Lambda(3) = ", L3, "   (= Lambda(1) : ", L1, ")");
xi = 4*Pi^3*L3;
print("F_4a : L(Xi_a,2) = L(Phi,3) = 4 pi^3 Lambda(3) = ", xi);
print("   xi/zeta(2)  = ", xi/zeta(2));
print("   xi/zeta(3)  = ", xi/zeta(3));
print("   xi/Pi^2     = ", xi/Pi^2);
print("   xi/Pi^3     = ", xi/Pi^3);
print("   lindep[xi,1,zeta(2),zeta(3)]      = ", lindep([xi,1,zeta(2),zeta(3)]));
print("   lindep[xi,zeta(2)]                = ", lindep([xi,zeta(2)]));
print("   lindep[xi,1]                      = ", lindep([xi,1]));
print("   lindep[xi,Pi^2,Pi^3,1]            = ", lindep([xi,Pi^2,Pi^3,1]));
print("   algdep(xi,4)                      = ", algdep(xi,4));
print("   algdep(xi/Pi^2,4)                 = ", algdep(xi/Pi^2,4));
print("   algdep(xi/Pi^3,4)                 = ", algdep(xi/Pi^3,4));
}
print("");
print("F_4b: the double pole sits at tau = i, i.e. ON the Fricke geodesic (0,i infty);");
print("      Phi_b(iy) ~ const/(y-1)^2, so int_0^inf y^{s-1} Phi_b(iy) dy DIVERGES.");
print("      |c(m)|^{1/m} -> exp(2 pi) exactly cancels the e^{-2 pi m} in the split integral.");
print("      partial sums of sum c(m)(2 pi m)^{-3} Gamma(3,2 pi m) at m = 100,200,300,400:");
{
my(cvb, s);
cvb = vector(LIM, m, polcoeff(Pb,m));
for(t=1,4, s = 100*t; print("        m<=", s, " : ", sum(m=1, s, cvb[m]*(2*Pi*m)^(-3)*incgam(3, 2*Pi*m))));
}
quit;
