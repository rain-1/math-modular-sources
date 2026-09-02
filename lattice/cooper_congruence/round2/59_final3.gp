\\ 59_final3.gp -- s18 identity with a CORRECTED Heegner-representative search.
\\ BUG in heeg.gp/heegrep: for a first column (p,r) it tests only ONE completion (q,s);
\\ the general matrix is [p, q+t p; r, s+t r], which shifts B by 2At.  Testing only t=0
\\ discards valid representatives and, in particular, the ones with the SMALLEST A -- which
\\ is why the s18 evaluations at the polar CM classes (content m, reduced disc -36) blew up.
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 150);
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
\\ minimal-A Heegner representative, all completions t
{ heegmin2(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3],d=b^2-4*a*c,A,B,B0,C,g,s,q,e,t,best=0);
  for(r=0,PB, for(p=-PB,PB,
    if(gcd(p,r)!=1, next);
    A = a*p^2+b*p*r+c*r^2;
    if(A<=0 || A%N!=0, next);
    if(best!=0 && A>=best[1], next);
    g = bezout(p,-r); if(g[3]!=1, next);
    s = g[1]; q = g[2];
    B0 = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
    \\ solve B0 + 2*A*t = beta (mod 2N)
    e = gcd(2*A, 2*N);
    if((beta-B0)%e != 0, next);
    t = lift(Mod((beta-B0)/e, 2*N/e) / Mod(2*A/e, 2*N/e));
    B = B0 + 2*A*t;
    C = (B^2-d)/(4*A);
    best = [A,B,C]));
  best;
}
initfser(3,1200);
{ fh(t) = if(abs(Fmod(3,t)) < 1e-40, fhatQ(t,18,FSER[3],#Vec(FSER[3][1])-3), fhatC(3,t)[2]); }
b18 = read("20_beta_s18.txt");
print("s18:  -1/(4 sqrt3) = ", -1/(4*sqrt(3)));
print("      -1/(2 sqrt3) = ", -1/(2*sqrt(3)));
print();
{
forstep(m=1,25,2,
  my(d=-36*m^2, bt=(18*m)%36, RF=redforms(d), T=0., nf=0, nadm=0, rep, om, al, Amx=0);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    nadm++;
    rep = heegmin2(RF[i],18,bt,60);
    if(rep==0, print("   NOREP ",RF[i]); next);
    nf++; om = omeg(rep); Amx = max(Amx,rep[1]);
    al = (-rep[2] + I*6*m)/(2*rep[1]);
    T += genchar(rep,-3)*fh(al)/om);
  print("  m=",m,"  cls=",nf,"/",nadm,"/",#RF,"  maxA=",Amx,"  beta=",b18[m]);
  print("      beta/T = ",b18[m]/T);
  print("      target = ",if(m%3, -1/(4*sqrt(3)), -1/(2*sqrt(3)))));
}
quit;
