\\ 58_s18fib2.gp -- is the Gamma_0(18)-Heegner fibre > 1 at NON-SQUARE d?  If so the trace
\\ over SL_2(Z)-classes is incomplete, which would explain the irrational values.
read("56_cdlib.gp");
default(realprecision, 50);
initfser(3,2200);
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
{ allreps2(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3],dd=b^2-4*a*c,L=List(),A,B,B0,C,g,s,q,e,t);
  for(r=0,PB, for(p=-PB,PB,
    if(gcd(p,r)!=1, next);
    A = a*p^2+b*p*r+c*r^2;
    if(A<=0 || A%N!=0, next);
    g = bezout(p,-r); if(g[3]!=1, next);
    s = g[1]; q = g[2];
    B0 = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
    e = gcd(2*A,2*N);
    if((beta-B0)%e != 0, next);
    t = lift(Mod((beta-B0)/e, 2*N/e) / Mod(2*A/e, 2*N/e));
    B = B0 + 2*A*t; C = (B^2-dd)/(4*A);
    listput(L,[A,B,C])));
  Vec(L);
}
{ dv(V,tol) = my(D=List()); for(i=1,#V, my(nw=1); for(j=1,#D, if(abs(V[i]-D[j])<tol*(1+abs(D[j])), nw=0;break)); if(nw, listput(D,V[i]))); Vec(D); }
CASES=[[5,30],[7,6],[13,6],[41,6],[9,18],[25,18]];
{
for(jj=1,#CASES,
  my(d=CASES[jj][1], bt=CASES[jj][2], D=-36*d, RF=redforms(D), tot=0, sumall=0.);
  print("=== d=",d,"  beta=",bt,"  disc=",D, if(issquare(d),"  (SQUARE)",""));
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    my(RS=allreps2(RF[i],18,bt,30), VL, DD, ch);
    if(#RS==0, print("   ",RF[i],"  NO REP"); next);
    VL = vector(#RS, j, fhx(3, (-RS[j][2]+I*sqrt(-D))/(2*RS[j][1])));
    DD = dv(VL, 1e-20);
    ch = genchar(RS[1],-3)*kronecker(-4,cont(RF[i]));
    tot += #DD;
    for(j=1,#DD, sumall += ch*DD[j]/omeg(RS[1]));
    print("   ",RF[i]," cont=",cont(RF[i]),"  chi*=",ch,"  reps=",#RS,"  DISTINCT Gamma_0(18)-classes = ",#DD));
  print("   total Gamma_0(18)-classes = ",tot,"   sum over ALL of them / sqrt3 = ",sumall/sqrt(3));
);
}
quit;
