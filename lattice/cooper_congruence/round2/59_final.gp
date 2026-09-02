\\ 59_final.gp -- high-precision verification of the fhat trace identities.
\\  s10: beta(m) = i*T(m)  (m odd),  = 2i*T(m)  (m even),  m<=30, 5 nmid m.
\\  s18: beta(m) = -(1/(4 sqrt3))*T_{chi_-3}(m)   (m odd, 3 nmid m)
\\       beta(m) = -(1/(2 sqrt3))*T_{chi_-3}(m)   (m odd, 3 | m),   m<=21.
\\ Heegner representatives are chosen with MINIMAL A (maximal Im) -- essential for s18.
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 150);
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
{ heegmin(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3],A,B,C,g,s,q,best=0);
  for(r=0,PB, for(p=-PB,PB,
    if(gcd(p,r)!=1, next);
    A = a*p^2+b*p*r+c*r^2;
    if(A<=0 || A%N!=0, next);
    if(best!=0 && A>=best[1], next);
    g = bezout(p,-r); if(g[3]!=1, next);
    s = g[1]; q = g[2];
    B = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
    if((B-beta)%(2*N)!=0, next);
    C = (B^2-(b^2-4*a*c))/(4*A);
    best = [A,B,C]));
  best;
}
initfser(2,900); initfser(3,900);
{ fh(k,t) = my(Fv=Fmod(k,t));
  if(abs(Fv) < 1e-40, return(fhatQ(t,LEV[k],FSER[k],#Vec(FSER[k][1])-3)));
  fhatC(k,t)[2];
}
b10 = read("20_beta_s10.txt");
b18 = read("20_beta_s18.txt");
print("======== s10:  beta(m) / (i T(m)),  T = sum chi_{-4}(Q) fhat(alpha_Q)/omega_Q ========");
{
for(m=1,30,
  if(m%5==0, print("  m=",m,"   anomalous (5|m)"); next);
  my(d=-4*m^2, bt=(6*m)%20, RF=redforms(d), T=0., nf=0, rep, om, al);
  for(i=1,#RF,
    rep = heegmin(RF[i],10,bt,60);
    if(rep==0, print("   NOREP ",RF[i]); next);
    nf++; om = omeg(rep);
    al = (-rep[2] + I*2*m)/(2*rep[1]);
    T += genchar(rep,-4)*fh(2,al)/om);
  print("  m=",m,"  cls=",nf,"/",#RF,"  beta=",b10[m],"   beta/(iT) = ", b10[m]/(I*T),
        "   predicted ",if(m%2,1,2)));
);
}
print();
print("======== s18:  beta(m) / T_{chi_-3}(m),  admissible forms (3 nmid content) ========");
print("   -1/(4 sqrt3) = ", -1/(4*sqrt(3)), "     -1/(2 sqrt3) = ", -1/(2*sqrt(3)));
{
forstep(m=1,21,2,
  my(d=-36*m^2, bt=(18*m)%36, RF=redforms(d), T=0., nf=0, nadm=0, rep, om, al);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    nadm++;
    rep = heegmin(RF[i],18,bt,60);
    if(rep==0, print("   NOREP ",RF[i]); next);
    nf++; om = omeg(rep);
    al = (-rep[2] + I*6*m)/(2*rep[1]);
    T += genchar(rep,-3)*fh(3,al)/om);
  print("  m=",m,"  cls=",nf,"/",nadm,"/",#RF,"  beta=",b18[m]);
  print("      T = ",T);
  print("      beta/T = ",b18[m]/T,"    predicted ",if(m%3, -1/(4*sqrt(3)), -1/(2*sqrt(3)))));
}
quit;
