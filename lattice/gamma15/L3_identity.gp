/* Independent check of the weight-3 companion identity found in task2 4(g):
      Re L(3,psi_4) = phi^5 * Im L(3,psi_4)
   for psi_4 the odd quartic character mod 5 with psi_4(2) = i.
   Also re-derives the two character sums behind the far-cusp period -11 zeta(2)/5. */
default(realprecision, 210);
ps4(n) = if(n%5==0, 0, if(n%5==1, 1, if(n%5==2, I, if(n%5==3, -I, -1))));
Lchi(ch, s) = 5^(-s)*sum(a=1,4, ch(a)*zetahurwitz(s, a/5));
ph5 = (11+5*sqrt(5))/2;
L3 = Lchi(ps4, 3); R3 = real(L3); I3 = imag(L3);
print("Re L(3,psi4) = ", R3);
print("Im L(3,psi4) = ", I3);
print("Re/Im        = ", R3/I3);
print("phi^5        = ", ph5);
print("|Re - phi^5 Im| = ", abs(R3 - ph5*I3));
print("lindep [ReL3, ImL3, sqrt5*ImL3] = ", lindep([R3, I3, sqrt(5)*I3], 60));
print("");
print("For contrast, at s = 2 (the fold-regularity ratio):");
L2 = Lchi(ps4,2); print("  Re L(2,psi4) - phi^5 Im L(2,psi4) = ", real(L2) - ph5*imag(L2), "  (nonzero: it IS -xi)");
print("  Phi_new at cusp 1/2 would need Re L3 = -phi^-5 Im L3:  -2 ImL3 - 2 phi^5 ReL3 = ", -2*I3 - 2*ph5*R3);
print("");
print("The two character sums behind the far-cusp period of Phi_D:");
re4(n) = real(ps4(n));
im4(n) = imag(ps4(n));
ww(n) = re4(n) - 2*im4(n);
print("  w(1..5)  = ", vector(5,j,ww(j)), "   sum_j j w(j) = ", sum(j=1,5, j*ww(j)), "  ->  ", -sum(j=1,5,j*ww(j))/5);
gg(n) = if(n%2==0, 1, -1/2);
print("  (w g)(1..10) = ", vector(10,j,ww(j)*gg(j)), "   sum_j j (wg)(j) = ", sum(j=1,10, j*ww(j)*gg(j)), "  ->  ", -sum(j=1,10,j*ww(j)*gg(j))/10);
print("  so pi^{(t1)}(Phi_D) = zeta(2)/5 and pi^{(t2)}(Phi_D) = -11 zeta(2)/5.");
quit;
