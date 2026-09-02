\\ 59_fibre10.gp -- TEST of hypothesis (ii): is  {Q : 10|A, B=beta (20), disc d}/Gamma_0(10)
\\ -> {disc d}/SL_2(Z)  a bijection?  fhat is Gamma_0(10)-invariant, so the number of
\\ DISTINCT fhat values over all Heegner representatives of one SL_2(Z)-class is a lower
\\ bound for the size of the fibre.
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 40);
N=10; initfser(2,700);
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
\\ all Heegner reps of Q in the (p,r) box
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
{ distinctvals(V, tol) = my(D=List());
  for(i=1,#V, my(new=1);
    for(j=1,#D, if(abs(V[i]-D[j])<tol, new=0; break));
    if(new, listput(D,V[i])));
  Vec(D);
}
{
MS=[11,13,14,16,12,17];
for(ii=1,#MS,
  my(m=MS[ii], d=-4*m^2, bt=(6*m)%20, RF=redforms(d));
  print("=== m=",m,"  d=",d,"  beta=",bt,"  #SL2 classes=",#RF);
  for(i=1,#RF,
    my(RS=allreps(RF[i],N,bt,26), VL, DV);
    if(#RS==0, print("   ",RF[i],"  no rep in box"); next);
    VL = vector(#RS, j, my(rp=RS[j]); fhatR(2, (-rp[2]+I*2*m)/(2*rp[1])));
    DV = distinctvals(VL, 1e-15);
    print("   ",RF[i]," cont=",cont(RF[i]),"  #reps in box=",#RS,"  #distinct fhat values=",#DV,
          if(#DV>1, concat("   *** FIBRE >1 ***  values: ",Str(DV)), "")));
);
}
quit;
