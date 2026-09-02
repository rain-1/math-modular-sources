\\ 59_s18fix.gp -- the CORRECTED s18 genus character.
\\ The discrepancy at m = 7,11,19,21,23 is exactly twice the content-p part of the trace,
\\ i.e. the classes of content g must be weighted by an extra chi_{-4}(g).  -36 is NOT
\\ fundamental: -36 = (-3)*12 = (-4)*9, and the character that works is the PRODUCT of the
\\ chi_{-3} genus character of the form with the chi_{-4} symbol of its content:
\\        chi*(Q) = chi_{-3}(Q) * kronecker(-4, cont(Q)).
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 120);
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
{ heegmin2(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3],d=b^2-4*a*c,A,B,B0,C,g,s,q,e,t,best=0);
  for(r=0,PB, for(p=-PB,PB,
    if(gcd(p,r)!=1, next);
    A = a*p^2+b*p*r+c*r^2;
    if(A<=0 || A%N!=0, next);
    if(best!=0 && A>=best[1], next);
    g = bezout(p,-r); if(g[3]!=1, next);
    s = g[1]; q = g[2];
    B0 = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
    e = gcd(2*A,2*N);
    if((beta-B0)%e != 0, next);
    t = lift(Mod((beta-B0)/e, 2*N/e) / Mod(2*A/e, 2*N/e));
    B = B0 + 2*A*t; C = (B^2-d)/(4*A);
    best = [A,B,C]));
  best;
}
initfser(3,1400);
{ fh(t) = if(abs(Fmod(3,t)) < 1e-40, fhatQ(t,18,FSER[3],#Vec(FSER[3][1])-3), fhatC(3,t)[2]); }
b18 = read("20_beta_s18.txt");
print("target lambda:  -1/(4 sqrt3) (3 nmid m),  -1/(2 sqrt3) (3 | m)");
{
forstep(m=1,29,2,
  my(d=-36*m^2, bt=(18*m)%36, RF=redforms(d), T=0., rep, om, al, tgt, nf=0);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    rep = heegmin2(RF[i],18,bt,60);
    if(rep==0, next);
    nf++; om = omeg(rep);
    al = (-rep[2] + I*6*m)/(2*rep[1]);
    T += genchar(rep,-3)*kronecker(-4,cont(RF[i]))*fh(al)/om);
  tgt = if(m%3, -1/(4*sqrt(3)), -1/(2*sqrt(3)));
  print("  m=",m," cls=",nf,"  beta/T - target = ", b18[m]/T - tgt, "     beta/T=",b18[m]/T));
}
quit;
