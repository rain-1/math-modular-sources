\\ 18_pzcontrol.gp -- CONTROL: does the same Maass-raised CM trace formula compute beta
\\ for Pasol-Zudilin's level-one magnetic forms?  F_4a = Delta/E_4^2 (double pole at rho,
\\ disc -3);  weight -2 weakly holomorphic form on SL_2(Z): f = E_4 E_6/Delta = q^{-1}-240-...
read("heeg.gp"); read("e2.gp");
default(realprecision, 60);
M = 260;
E4s = 1 + 240*sum(n=1,M,sigma(n,3)*'q^n) + O('q^(M+1));
E6s = 1 - 504*sum(n=1,M,sigma(n,5)*'q^n) + O('q^(M+1));
Dl  = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
F4a = Dl/E4s^2;
F4b = E4s*Dl/E6s^2;
fw  = E4s*E6s/Dl;
Dfw = 'q*deriv(fw,'q);
{ cvec(S) = vector(M,m,polcoeff(S,m)); }
{ betav(c) = my(cp=vector(M,m,c[m]/m), b=vector(M));
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*cp[n/d]); b[n]=s); b; }
B4a = betav(cvec(F4a));
B4b = betav(cvec(F4b));
print("F_4a: c(1..6) = ", vector(6,i,polcoeff(F4a,i)));
print("beta_4a(1..8) = ", vector(8,i,B4a[i]));
print("a(m^2) = m beta(m) : ", vector(6,i,i*B4a[i]));
print("f = E4E6/Delta : ", vector(8,i,polcoeff(fw,i-2)));
\\ SL_2(Z) reduction and evaluation of fhat
{ evs(f,q,TT) = my(v=valuation(f,'q), s=0.); for(n=v,TT, s += polcoeff(f,n)*q^n); s; }
{ fhat1(t) = my(R=sl2red(t), z, q, y);
  z = R[1]; y = imag(z); q = exp(2*Pi*I*z);
  evs(Dfw,q,M-3) + evs(fw,q,M-3)/(2*Pi*y);
}
{ tr(m, tw) = my(d=-3*m^2, RF=redforms(d), t=0., ch, om, al);
  for(i=1,#RF,
    ch = if(tw, genchar(RF[i],-3), 1);
    om = omeg(RF[i]);
    al = (-RF[i][2] + I*m*sqrt(3))/(2*RF[i][1]);
    t += ch*fhat1(al)/om);
  t;
}
{
for(m=1,14,
  my(A=tr(m,1), B=tr(m,0), a, b);
  a = real(I*sqrt(3)*A/(I*sqrt(3)));
  b = real(I*sqrt(3)*B/(I*sqrt(3)));
  print("m=",m,"  Tr(chi)=", A/I, "   Tr(nochi)=", B/I, "   beta_4a=", B4a[m], "   Tr(chi)/(192 I) - beta = ", A/(192*I) - B4a[m]);
);
}
quit;
