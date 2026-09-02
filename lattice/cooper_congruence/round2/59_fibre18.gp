\\ 59_fibre18.gp -- for s18 at the four deviant m: is the Gamma_0(18)-Heegner fibre over an
\\ SL_2(Z)-class bigger than 1?  fhat is Gamma_0(18)-invariant, so distinct fhat values over
\\ the reps of one SL_2(Z)-class = distinct Gamma_0(18)-classes.
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N=18; initfser(3,900);
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
{ allreps(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3],L=List(),A,B,C,g,s,q);
  for(r=0,PB, for(p=-PB,PB,
    if(gcd(p,r)!=1, next);
    A = a*p^2+b*p*r+c*r^2;
    if(A<=0 || A%N!=0, next);
    g = bezout(p,-r); if(g[3]!=1, next);
    s = g[1]; q = g[2];
    B = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
    if((B-beta)%(2*N)!=0, next);
    C = (B^2-(b^2-4*a*c))/(4*A);
    listput(L,[A,B,C])));
  Vec(L);
}
{ fh(t) = if(abs(Fmod(3,t)) < 1e-30, fhatQ(t,18,FSER[3],#Vec(FSER[3][1])-3), fhatC(3,t)[2]); }
{ dv(V,tol) = my(D=List()); for(i=1,#V, my(nw=1); for(j=1,#D, if(abs(V[i]-D[j])<tol*(1+abs(D[j])), nw=0;break)); if(nw, listput(D,V[i]))); Vec(D); }
MS=[5,7,11,13];
{
for(jj=1,#MS,
  my(m=MS[jj], d=-36*m^2, bt=(18*m)%36, RF=redforms(d), TT=0.);
  print("=== m=",m,"  d=",d,"  beta=",bt);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    my(RS=allreps(RF[i],N,bt,30), VL, D, Amin=10^9);
    if(#RS==0, print("   ",RF[i]," NO REP IN BOX"); next);
    for(j=1,#RS, Amin=min(Amin,RS[j][1]));
    VL = vector(#RS, j, fh((-RS[j][2]+I*6*m)/(2*RS[j][1])));
    D = dv(VL, 1e-25);
    if(#D>1,
      print("   ",RF[i]," cont=",cont(RF[i])," minA=",Amin," reps=",#RS,"  *** ",#D," DISTINCT CLASSES ***  ",D),
      print("   ",RF[i]," cont=",cont(RF[i])," minA=",Amin," reps=",#RS,"  1 class, fhat=",D[1])));
);
}
quit;
