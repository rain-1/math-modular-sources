default(parisizemax,6000000000);
default(realprecision,80);
M=600; default(seriesprecision,M+6);
et(k)=eta(q^k+O(q^(M+6)));
F=et(2)^7*et(3)^7/(et(1)^5*et(6)^5);
t=q*(et(1)*et(6)/(et(2)*et(3)))^12;
f=sqrt(F); Psi = f^3*t/4;
ps=vector(M,m,polcoeff(Psi,m));
xi=0.1001874492293394061677586821330611216287777205008122117162193987349722087120577580850958409788906881;
qs = exp(-2*Pi/sqrt(6));
Th = sum(m=1,M, ps[m]*qs^m/m^2);
print("Theta(q_*)       = ",Th);
print("Theta(q_*) - xi  = ",Th-xi);
print("Theta(q_*) / xi  = ",Th/xi);
\\ folded Mellin:  Lambda(Psi,s) = int_{1/sqrt6}^inf Psi(iy)[y^{s-1} + 6^{3/2-s} y^{2-s}] dy
Pv(y) = sum(m=1,M, ps[m]*exp(-2*Pi*m*y));
Lam2 = intnum(y=1/sqrt(6), [oo,2*Pi], Pv(y)*(y + 6^(-1/2)));
L2 = 4*Pi^2*Lam2;
print("L(Psi,2) = ",L2);
print("L(Psi,2) - xi = ",L2-xi);
print("L(Psi,2)/xi   = ",L2/xi);
\\ L(Psi,1) and the forced annihilation checks
Lam1 = intnum(y=1/sqrt(6), [oo,2*Pi], Pv(y)*(1 + 6*y));
print("L(Psi,1) = ",2*Pi*Lam1);
Lam3 = intnum(y=1/sqrt(6), [oo,2*Pi], Pv(y)*(y^2 + 6^(-3/2)/y));
print("L(Psi,3) = ",(2*Pi)^3/2*Lam3);
\q
