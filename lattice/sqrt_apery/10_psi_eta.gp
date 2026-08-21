default(parisizemax,6000000000);
M=200; default(seriesprecision,M+8);
et(k)=eta(q^k+O(q^(M+8)));
F=et(2)^7*et(3)^7/(et(1)^5*et(6)^5);
t=q*(et(1)*et(6)/(et(2)*et(3)))^12;
f=sqrt(F); Psi=f^3*t/4;
\\ claim:  f  = (eta2 eta3)^{7/2}/(eta1 eta6)^{5/2}     [q-order 0]
\\ claim:  Psi= (1/4)(eta1 eta6)^{9/2}/(eta2 eta3)^{3/2} [q-order 1]
fclaim = sqrt((et(2)*et(3))^7/(et(1)*et(6))^5);
Pclaim = q*sqrt((et(1)*et(6))^9/((et(2)*et(3))^3))/4;
print("f  == (eta2 eta3)^{7/2}/(eta1 eta6)^{5/2} ? ", f-fclaim==0);
print("Psi== (1/4)(eta1 eta6)^{9/2}/(eta2 eta3)^{3/2} ? ", Psi-Pclaim==0);
print("Psi coeffs*4^m (integers) m=1..14: ",vector(14,m,4^m*polcoeff(Psi,m)));
print("v_2(den psi(m)) m=1..60: ",vector(60,m,valuation(denominator(polcoeff(Psi,m)),2)));
\q
