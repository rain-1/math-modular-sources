\\ 19_pzcontrol2.gp -- level-one control, both PZ forms:
\\   F_4a = Delta/E_4^2 (double pole at rho, D_0=-3),  F_4b = E_4 Delta/E_6^2 (at i, D_0=-4).
\\ f = E_4E_6/Delta = q^{-1}-240-... is THE weight -2 weakly holomorphic form on SL_2(Z)
\\ with pole order 1.  fhat = Df + f/(2 pi y).
read("heeg.gp"); read("e2.gp");
default(realprecision, 60);
M = 300;
E4s = 1 + 240*sum(n=1,M,sigma(n,3)*'q^n) + O('q^(M+1));
E6s = 1 - 504*sum(n=1,M,sigma(n,5)*'q^n) + O('q^(M+1));
Dl  = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
F4a = Dl/E4s^2;
F4b = E4s*Dl/E6s^2;
F6  = E6s*Dl/E4s^3;
fw  = E4s*E6s/Dl;
Dfw = 'q*deriv(fw,'q);
{ betav(S) = my(c=vector(M,m,polcoeff(S,m)), cp, b=vector(M));
  cp = vector(M,m,c[m]/m);
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*cp[n/d]); b[n]=s); b; }
B4a = betav(F4a); B4b = betav(F4b);
{ evs(f,q,TT) = my(v=valuation(f,'q), s=0.); for(n=v,TT, s += polcoeff(f,n)*q^n); s; }
{ fhat1(t) = my(R=sl2red(t), z, q, y);
  z = R[1]; y = imag(z); q = exp(2*Pi*I*z);
  evs(Dfw,q,M-3) + evs(fw,q,M-3)/(2*Pi*y);
}
{ tr(m, D0) = my(d=D0*m^2, RF=redforms(d), t=0.);
  for(i=1,#RF, t += genchar(RF[i],D0)*fhat1((-RF[i][2]+I*m*sqrt(-D0))/(2*RF[i][1]))/omeg(RF[i]));
  t;
}
print("F_4a (D_0 = -3):   beta(m)  vs  I*Tr/192");
{ for(m=1,12, my(T=tr(m,-3)); print("  m=",m,"  beta=",B4a[m],"   I*Tr/192 - beta = ", I*T/192 - B4a[m])); }
print("");
print("F_4b (D_0 = -4):   beta(m)  vs  lam*Tr,  lam fitted from m=1");
{ my(T1=tr(1,-4), lam=B4b[1]/T1);
  print("  lam = ", lam, "   (2I/lam = ", 2*I/lam, ")");
  for(m=1,12, my(T=tr(m,-4)); print("  m=",m,"  beta=",B4b[m],"   lam*Tr - beta = ", lam*T - B4b[m])); }
print("");
print("divisibility check, level one:  v_2(beta_4a(2))=",valuation(B4a[2],2),"  beta_4a(2)=",B4a[2]);
print("  m^2 | beta_4a(m)?  ", vector(8,m,if(m==1,1,B4a[m]%(m^2)==0)));
quit;
