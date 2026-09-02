\\ 51_master.gp -- TASK 3: master verification n^2 | beta(n) for n <= M=3000, all rows;
\\ sharpness v_p(gamma(n)); zeros of gamma; refined supercongruences of COOPER_CONGRUENCE.md sec.4.
read("50_lib.gp");
M = 12000;
PL = primes([2,60]);
gettime();
GG = vector(3); BB = vector(3); CP = vector(3);
{
for(k=1,3,
  my(cp, b, g, bad=0);
  cp = CPvec(k,M);
  b  = Bvec(k,cp);
  g  = vector(M);
  for(n=1,M, if(b[n]%(n^2)!=0, bad=n; break, g[n]=b[n]/n^2));
  CP[k]=cp; BB[k]=b; GG[k]=g;
  print("row ", NAM[k], "  n^2 | beta(n) for n<=", M, ": ", if(bad==0,"YES",concat("FAIL at ",bad)), "   [", gettime(), " ms]");
  ;
);
}
print();
print("== zeros of gamma(n), n <= ", M, " ==");
{
for(k=1,3,
  my(Z=select(n->GG[k][n]==0, vector(M,i,i)));
  print("row ",NAM[k],": #zeros=",#Z,"   are they exactly the multiples of ",if(k==1,7,if(k==2,5,0)),"? ",
        if(k==3, #Z==0, Z == vector(#Z,i,i*if(k==1,7,5))));
);
}
print();
print("== sharpness: min over n<=M of v_p(gamma(n)) (= min of v_p(beta(n)) - 2 v_p(n)) ==");
print("   A_p = min over all n with gamma(n)!=0 ;  B_p = min over n with p|n, gamma(n)!=0");
{
for(k=1,3,
  print("row ",NAM[k]);
  for(i=1,#PL,
    my(p=PL[i], a=oo, bb=oo, na=0, nb=0, v);
    for(n=1,M,
      if(GG[k][n]==0, next);
      v = valuation(GG[k][n],p);
      if(v<a, a=v; na=n);
      if(n%p==0 && v<bb, bb=v; nb=n));
    print("  p=",p,"  A_p=",a," (at n=",na,")   B_p=",bb," (at n=",nb,")"));
  print();
);
}
print("== refined supercongruence 1: c'(p^k) = psi(p) c'(p^(k-1)) mod p^(2k),  p^k <= ",M," ==");
{
for(k=1,3,
  my(fails=List(), tot=0);
  forprime(p=2,M,
    my(e=1, ps=psival(k,p));
    while(p^e<=M,
      tot++;
      if((CP[k][p^e] - ps*if(e==1,1,CP[k][p^(e-1)])) % p^(2*e) != 0,
         listput(fails,[p,e]));
      e++));
  print("row ",NAM[k],": tested ",tot," (p,k) pairs; failures: ",if(#fails==0,"NONE",Vec(fails)));
);
}
print();
print("SKIP1b"); print("== refined supercongruence 1b: exact modulus reached, p^k <= ",M,", k>=1 ==");
{
for(k=1,3,
  print("row ",NAM[k],":");
  forprime(p=2,53,
    my(e=1, ps=psival(k,p), s="");
    while(p^e<=M,
      my(D = CP[k][p^e] - ps*if(e==1,1,CP[k][p^(e-1)]), v);
      v = if(D==0, oo, valuation(D,p));
      s = concat(s, concat(concat(concat("  k=",e)," v_p(diff)="),concat(Str(v),concat("/need ",Str(2*e)))));
      e++);
    print("   p=",p,s)));
}
print();
print("== refined supercongruence 2 (two primes): c'(pq)-psi(p)c'(q)-psi(q)c'(p)+psi(pq) = 0 mod p^2 q^2 ==");
{
for(k=1,3,
  my(fails=List(), tot=0);
  forprime(p=2,M,
    forprime(q=p+1,M\p,
      my(D);
      tot++;
      D = CP[k][p*q] - psival(k,p)*CP[k][q] - psival(k,q)*CP[k][p] + psival(k,p)*psival(k,q);
      if(D % (p^2*q^2) != 0, listput(fails,[p,q]))));
  print("row ",NAM[k],": tested ",tot," pairs p<q, pq<=",M,"; failures: ",if(#fails==0,"NONE",Vec(fails)));
);
}
print();
print("== generalised: c'(n) = sum_{d|n} psi(n/d) d^2 gamma(d)  (Conj 4.1) for n <= ",M," ==");
{
for(k=1,3,
  my(bad=0);
  for(n=1,M,
    my(s=0);
    fordiv(n,d, s += psin(k,n/d)*d^2*GG[k][d]);
    if(s != CP[k][n], bad=n; break));
  print("row ",NAM[k],": ",if(bad==0,"HOLDS for all n<=M",concat("FAILS at n=",bad)));
);
}
print();
print("total time ", gettime(), " ms");
quit;
