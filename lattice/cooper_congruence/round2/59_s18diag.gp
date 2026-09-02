\\ 59_s18diag.gp -- locate the s18 discrepancy at m = 7, 11, 19, 21, 23.
\\ (Successes: m = 1,3,5,9,13,15,17,25;  failures: 7,11,19,21,23 -- exactly the m with a
\\  prime factor p = 3 mod 4, p != 3, i.e. a prime INERT in Q(i) other than the conductor 3.)
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 90);
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
    e = gcd(2*A, 2*N);
    if((beta-B0)%e != 0, next);
    t = lift(Mod((beta-B0)/e, 2*N/e) / Mod(2*A/e, 2*N/e));
    B = B0 + 2*A*t; C = (B^2-d)/(4*A);
    best = [A,B,C]));
  best;
}
initfser(3,1200);
{ fh(t) = if(abs(Fmod(3,t)) < 1e-40, fhatQ(t,18,FSER[3],#Vec(FSER[3][1])-3), fhatC(3,t)[2]); }
b18 = read("20_beta_s18.txt");
MS=[7,11,19,21,23,5,13];
{
for(jj=1,#MS,
  my(m=MS[jj], d=-36*m^2, bt=(18*m)%36, RF=redforms(d), T=0., Tgt, byc=Map(), rep, om, al, tm, cg);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    rep = heegmin2(RF[i],18,bt,60);
    if(rep==0, next);
    om = omeg(rep);
    al = (-rep[2] + I*6*m)/(2*rep[1]);
    tm = genchar(rep,-3)*fh(al)/om;
    T += tm;
    cg = cont(RF[i]);
    mapput(byc, cg, if(mapisdefined(byc,cg), mapget(byc,cg), 0.) + tm));
  Tgt = b18[m]/if(m%3, -1/(4*sqrt(3)), -1/(2*sqrt(3)));
  print("m=",m,"   T=",T,"   T_target=",Tgt);
  print("    (T - T_target)/sqrt(3) = ", (T-Tgt)/sqrt(3));
  print("    partial sums by content:");
  my(K=Vec(Mat(byc)[,1]));
  for(i=1,#K, print("        cont=",K[i],"   sum=",mapget(byc,K[i]),"   /sqrt3=",mapget(byc,K[i])/sqrt(3)));
);
}
quit;
