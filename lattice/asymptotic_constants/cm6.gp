default(parisize,"14G");
default(realprecision,200);
M = 950;
et(k)=prod(m=1,M\k+1,1-'q^(k*m)+O('q^(M+2)));
D(s)='q*deriv(s,'q);
ev(s,z)=subst(truncate(s),'q,z);
E8(k)=1+480*sum(m=1,M\k+1,sigma(m,7)*'q^(k*m))+O('q^(M+2));
h = et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/'q;
v = 1/h;
X = 'q*et(4)^2*et(12)^2/(et(1)^2*et(3)^2);
T = X/(16*X^2+2*X+1);
PH = (E8(1)-572*E8(2)+11583*E8(3)-36608*E8(4)+46332*E8(6)-20736*E8(12))/480;
F = PH/D(T);
qc = exp(-Pi/sqrt(3));
Dv = ev(D(v),qc); Fv = ev(F,qc); D2T = ev(D(D(T)),qc); Dx = ev(D(X),qc); vc = ev(v,qc);
print("v(tau_c)      = ", vc, "   1/(2 sqrt3) = ", 1/(2*sqrt(3)));
print("Dv(tau_c)     = ", Dv);
A1 = Fv/Dv^3; A2 = D2T/Dv^2; A3 = Dx/Dv;
print("F/(Dv)^3      = ", A1);
print("  algdep 2: ", algdep(A1,2), "  algdep 4: ", algdep(A1,4));
print("D2t/(Dv)^2    = ", A2);
print("  algdep 2: ", algdep(A2,2), "  algdep 4: ", algdep(A2,4));
print("Dx/Dv         = ", A3);
print("  algdep 2: ", algdep(A3,2), "  algdep 6: ", algdep(A3,6));
\\ rational fit of x_{zeta7} in v
rf(GG,U,d1,d2,NC)={
  my(cols=List());
  for(j=0,d1, listput(cols, Vec(truncate(U^j+O('q^(NC+1))),-(NC+1))));
  for(j=0,d2, listput(cols, Vec(truncate(-GG*U^j+O('q^(NC+1))),-(NC+1))));
  Mat(vector(#cols,i,cols[i]~));
}
{for(d=1,6, my(k=matker(rf(X,v,d,d,6*d+80))); if(#k>0, print("x_{z7} rational in v, degree ",d,": ",k~); break));}
{for(d=1,10, my(k=matker(rf(v^3*F/D(v)^3,v,d,d,6*d+120))); if(#k>0, print("v^3 F/(Dv)^3 rational in v, degree ",d,": ",k~); break));}
Kf = sqrt((1/10)*ev(D(F),qc)^2/(2*Pi*(-D2T)));
W2=gamma(1/3)^6/Pi^4;
kap = Kf*Pi^(3/2)/W2^2;
print("K                = ", Kf);
print("kappa = K pi^{3/2}/W2^2 = ", kap);
{foreach([1,2,3,4,6,8,12],dd, print("  algdep(kappa,",dd,") = ", algdep(kap,dd)));}
{foreach([2,3,6],e, print("  algdep(kappa^",e,",2) = ", algdep(kap^e,2), "   algdep(kappa^",e,",4) = ", algdep(kap^e,4)));}
print("  kappa^3 = ", kap^3, "  kappa^6 = ", kap^6);
print("  kappa^6 * 2^68 / 30^3 = ", kap^6*2^68/27000);
quit;
