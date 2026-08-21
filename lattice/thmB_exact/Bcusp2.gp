read("common.gp");
default(realprecision, 60);
MM = 900000;
cv = vector(MM, m, cm("B",m)*1.0/m^2);
Th(q0,pr) = my(M=min(MM,ceil(pr*log(10)/(-log(abs(q0))))), s=0.); forstep(m=M,1,-1, s=(s+cv[m])*q0); s;
G(Y) = my(sig=I*Y, tau=sig/(6*sig+1), q0=exp(2*Pi*I*tau)); Th(q0,55)*(6*sig+1);
{
Y1=18.; Y2=23.;
g1=G(Y1); g2=G(Y2);
c1=(g1-g2)/(I*Y1-I*Y2); c0=g1-c1*I*Y1;
xiB = c1/6;
print("c1 = ",c1);
print("c0 = ",c0);
print("xi_B = c1/6 = ", xiB);
print("L(Phi_B,2)  = ", Ltarget("B"));
print("Re xi_B - L = ", real(xiB)-Ltarget("B"));
print("Im xi_B     = ", imag(xiB));
print("2*Pi^2/(27*sqrt(3)) = ", 2*Pi^2/(27*sqrt(3)));
print("diff = ", imag(xiB)-2*Pi^2/(27*sqrt(3)));
print("lindep(Im, Pi^2/sqrt(3)) = ", lindep([imag(xiB), Pi^2/sqrt(3)]));
}
quit;
