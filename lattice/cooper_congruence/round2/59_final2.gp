\\ 59_final2.gp -- (a) the s10 identity to m<=30; (b) the s18 cases that were only
\\ numerically approximate in 59_final.log, redone with a larger search box (PB=200,
\\ so that the MINIMAL A -- hence the maximal Im -- is really found) and 250 digits.
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 250);
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
initfser(2,1200); initfser(3,1200);
{ fh(k,t) = if(abs(Fmod(k,t)) < 1e-60, fhatQ(t,LEV[k],FSER[k],#Vec(FSER[k][1])-3), fhatC(k,t)[2]); }
b10 = read("20_beta_s10.txt");
b18 = read("20_beta_s18.txt");
print("======== s10: beta(m)/(i T(m)) ; predicted 1 (m odd), 2 (m even) ========");
{
for(m=1,30,
  if(m%5==0, print("  m=",m,"   anomalous (5|m)"); next);
  my(d=-4*m^2, bt=(6*m)%20, RF=redforms(d), T=0., nf=0, rep, om, al, Amx=0);
  for(i=1,#RF,
    rep = heegmin(RF[i],10,bt,100);
    if(rep==0, print("   NOREP ",RF[i]); next);
    nf++; om = omeg(rep); Amx=max(Amx,rep[1]);
    al = (-rep[2] + I*2*m)/(2*rep[1]);
    T += genchar(rep,-4)*fh(2,al)/om);
  print("  m=",m,"  cls=",nf,"/",#RF,"  maxA=",Amx,"  pred=",if(m%2,1,2),"   beta/(iT) = ", b10[m]/(I*T)));
}
print();
print("======== s18: the four cases m = 7, 11, 19, 21 redone with PB=200, 250 digits ========");
MS=[7,11,19,21];
{
for(j=1,#MS,
  my(m=MS[j], d=-36*m^2, bt=(18*m)%36, RF=redforms(d), T=0., nf=0, nadm=0, rep, om, al, Amx=0);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    nadm++;
    rep = heegmin(RF[i],18,bt,200);
    if(rep==0, print("   NOREP ",RF[i]); next);
    nf++; om = omeg(rep); Amx=max(Amx,rep[1]);
    al = (-rep[2] + I*6*m)/(2*rep[1]);
    T += genchar(rep,-3)*fh(3,al)/om);
  print("  m=",m,"  cls=",nf,"/",nadm,"  maxA=",Amx,"  beta=",b18[m]);
  print("      beta/T   = ",b18[m]/T);
  print("      predicted= ",if(m%3, -1/(4*sqrt(3)), -1/(2*sqrt(3)))));
}
quit;
