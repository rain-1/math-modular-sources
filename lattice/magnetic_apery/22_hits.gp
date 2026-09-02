default(parisizemax,12000000000);
read("lib.gp");
NQ = 1005; MC = 1000;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
print("### the four magnetic hits on Gamma_0(6): exact data");
\\ [name, Q_u, alpha]
NAMES = ["FOLD: Q=(72u^2-1)^4", "CUSP: Q=(72u^2+17u+1)^3", "D8: Q=((72u^2+16u+1)(72u^2+17u+1))^2", "D12: Q=((72u^2+18u+1)(72u^2+17u+1))^2"];
QQ = [(72*x^2-1)^4, (72*x^2+17*x+1)^3, ((72*x^2+16*x+1)*(72*x^2+17*x+1))^2, ((72*x^2+18*x+1)*(72*x^2+17*x+1))^2];
AA = [[0,-2,-68,-864,-4896,-10368,0], [-5,-84,0,6048,25920], [-1,-35,-378,0,27216,181440,373248], [-1,-37,-544,-3972,-14400,-20736,0]];
{
for(t=1, 4,
  my(Qp, Nn, rho, Phi, Xi, bad, gg, dd);
  Qp = QQ[t];
  Nn = sum(k=0, #AA[t]-1, AA[t][k+1]*x^(k+1));
  gg = gcd(Nn, Qp);
  print("");
  print("== ", NAMES[t]);
  print("   numerator N(u) = ", Nn);
  print("   N factored     = ", factor(Nn));
  print("   gcd(N,Q)       = ", gg);
  print("   reduced rho    = ", Nn/gg, "  /  ", Qp/gg);
  rho = subst(Nn, x, us)/subst(Qp, x, us);
  Phi = F2*rho;
  Xi = Dinv(Phi, MC);
  print("   c(1..12) = ", vector(12,j,polcoeff(Phi,j)));
  print("   integrality of Phi: first failure m = ", intfail(Phi, MC));
  print("   MAGNETIC (m | c(m)): first failure m = ", magfail(Phi, MC), "   [0 = holds for all m <= ", MC, "]");
  print("   c'(1..12) = ", vector(12,j,polcoeff(Phi,j)/j));
  print("   |c(m)|^(1/m) at m=", MC, " : ", abs(polcoeff(Phi,MC)*1.0)^(1/MC));
);
}
print("");
print("reference growth rates: e^(2 pi Im tau_0) for Im = 1/sqrt6, sqrt2/6, 1/(4 sqrt3):");
print("  fold  1/sqrt6 = ", 1/sqrt(6), " -> ", exp(2*Pi/sqrt(6)));
print("  D-8   sqrt2/6 = ", sqrt(2)/6, " -> ", exp(2*Pi*sqrt(2)/6));
print("  D-12  1/(4sqrt3) = ", 1/(4*sqrt(3)), " -> ", exp(2*Pi/(4*sqrt(3))));
quit;
