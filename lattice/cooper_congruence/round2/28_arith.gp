\\ 28_arith.gp -- what is gamma(p) mod p ?  Sweep of arithmetic candidates.
read("lib.gp");
G  = [read("22_gamma_s7.txt"), read("22_gamma_s10.txt"), read("22_gamma_s18.txt")];
N  = #G[1];
PMAX = 200;
PL = select(x->isprime(x)&&x>7, vector(PMAX,i,i));
print("primes used: ", #PL, " up to ", PMAX);

\\ target residues
{ tgt(k) = vector(#PL, i, lift(Mod(G[k][PL[i]], PL[i]))); }
T = [tgt(1),tgt(2),tgt(3)];
print();
print("target gamma(p) mod p:");
for(k=1,3, print("  ",NAM[k],": ",T[k]));

\\ compare a candidate vector (indexed like PL) with the target, +/- and inverse
{ chk(name, C) =
  for(k=1,3,
    my(m1=0,m2=0,m3=0,m4=0,tot=0);
    for(i=1,#PL, my(p=PL[i]);
      if(C[i]===0 && 0, );
      tot++;
      if((T[k][i]-C[i])%p==0, m1++);
      if((T[k][i]+C[i])%p==0, m2++);
      if(C[i]%p!=0 && (T[k][i]*C[i]-1)%p==0, m3++);
      if(C[i]%p!=0 && (T[k][i]*C[i]+1)%p==0, m4++);
    );
    if(m1>4 || m2>4 || m3>4 || m4>4,
      print("  HIT  ",NAM[k]," vs ",name,"  = :",m1,"/",tot,"  = -:",m2,"/",tot,"  inv:",m3,"/",tot,"  -inv:",m4,"/",tot));
  );
}

print();
print("=== candidate sweep (a 'HIT' line is printed only if >4 of ",#PL," primes match) ===");

\\ (1) modular forms: trace forms of new spaces
{
for(k2=1,3, my(wt=[2,4,6][k2]);
  for(M=1,80,
    my(mf, tf, cf);
    mf = mfinit([M,wt],0);
    if(mfdim(mf)==0, next);
    tf = mftraceform([M,wt],0);
    cf = mfcoefs(tf, PMAX);
    chk(concat(concat(concat("Tr S_",wt),concat("^new(",M)),")"), vector(#PL,i,cf[PL[i]+1]));
  );
);
}

\\ (2) eta-quotient / CM candidates
{
my(M=PMAX, E1, E2, E3, cf);
E1 = Ed(1,M+2); E2 = Ed(2,M+2); E3 = Ed(3,M+2);
cf = Vec('q*E3^8 + O('q^(M+2)), -(M+1));   \\ eta(3t)^8 : weight 4 level 9
chk("eta(3t)^8", vector(#PL,i,cf[PL[i]+1]));
cf = Vec('q*(E1*Ed(7,M+2))^4 + O('q^(M+2)), -(M+1));
chk("eta(t)^4eta(7t)^4", vector(#PL,i,cf[PL[i]+1]));
cf = Vec('q*(E1*E2)^12/1 + O('q^(M+2)), -(M+1));
chk("eta(t)^12eta(2t)^12", vector(#PL,i,cf[PL[i]+1]));
cf = Vec('q*E1^24 + O('q^(M+2)), -(M+1));
chk("tau (Delta)", vector(#PL,i,cf[PL[i]+1]));
}

\\ (3) elementary arithmetic candidates
{
chk("constant 1",           vector(#PL,i,1));
chk("kronecker(-3,p)",      vector(#PL,i,kronecker(-3,PL[i])));
chk("kronecker(-4,p)",      vector(#PL,i,kronecker(-4,PL[i])));
chk("kronecker(-7,p)",      vector(#PL,i,kronecker(-7,PL[i])));
chk("kronecker(5,p)",       vector(#PL,i,kronecker(5,PL[i])));
chk("kronecker(-8,p)",      vector(#PL,i,kronecker(-8,PL[i])));
chk("B_{p-3} mod p",        vector(#PL,i,my(p=PL[i],b=bernfrac(p-3)); lift(Mod(numerator(b),p)/Mod(denominator(b),p))));
chk("B_{(p-1)/2} mod p",    vector(#PL,i,my(p=PL[i],b=bernfrac((p-1)\2)); if(denominator(b)%p==0,0,lift(Mod(numerator(b),p)/Mod(denominator(b),p)))));
chk("Fermat q_p(2)",        vector(#PL,i,my(p=PL[i]); lift(Mod((2^(p-1)-1)/p,p))));
chk("Fermat q_p(3)",        vector(#PL,i,my(p=PL[i]); lift(Mod((3^(p-1)-1)/p,p))));
chk("Fermat q_p(5)",        vector(#PL,i,my(p=PL[i]); lift(Mod((5^(p-1)-1)/p,p))));
chk("Fermat q_p(7)",        vector(#PL,i,my(p=PL[i]); lift(Mod((7^(p-1)-1)/p,p))));
chk("central binom",        vector(#PL,i,my(p=PL[i]); lift(Mod(binomial(p-1,(p-1)\2),p))));
chk("h(-p) or 0",           vector(#PL,i,my(p=PL[i]); if(p%4==3, qfbclassno(-p), 0)));
chk("h(-4p)",               vector(#PL,i,qfbclassno(-4*PL[i])));
chk("Wilson quotient",      vector(#PL,i,my(p=PL[i]); lift(Mod(((p-1)!+1)/p,p))));
chk("p mod 7",              vector(#PL,i,PL[i]%7));
chk("p mod 5",              vector(#PL,i,PL[i]%5));
chk("p mod 12",             vector(#PL,i,PL[i]%12));
chk("p mod 9",              vector(#PL,i,PL[i]%9));
chk("harmonic H_{(p-1)/2}", vector(#PL,i,my(p=PL[i],s=Mod(0,p)); for(j=1,(p-1)\2, s+=Mod(j,p)^-1); lift(s)));
chk("sum 1/j^2, j<=(p-1)/2",vector(#PL,i,my(p=PL[i],s=Mod(0,p)); for(j=1,(p-1)\2, s+=Mod(j,p)^-2); lift(s)));
}

\\ (3b) individual rational eigenforms in S_k^new(Gamma_0(M))
{
for(k2=1,3, my(wt=[2,4,6][k2]);
  for(M=1,60,
    my(mf, eb);
    mf = mfinit([M,wt],0);
    if(mfdim(mf)==0, next);
    eb = mfeigenbasis(mf);
    for(j=1,#eb,
      my(cf = mfcoefs(eb[j], PMAX), rat=1);
      for(i=1,#PL, if(type(cf[PL[i]+1])!="t_INT", rat=0; break));
      if(rat, chk(concat(concat(concat(concat("newform wt",wt),concat(" lev",M)),"#"),j), vector(#PL,i,cf[PL[i]+1])));
    );
  );
);
}

\\ (4) the x-side: a_{p-1} and its refinement
{
print();
print("=== x-side: eta = sum a_j x^j dx ; is gamma(p) = (a_{p-1}-psi(p))/p mod p ? ===");
for(k=1,3,
  my(R=ROWS[k], B=R[3], C=R[4], NM=PMAX+4, AB, A, Fx, lx, Px, aa, cf, ok1=0, ok2=0, tot=0);
  AB = genrow(k,NM);
  A  = AB[1];
  Fx = sum(n=0,NM,A[n+1]*'x^n) + O('x^(NM+1));
  lx = sum(n=0,NM,A[n+1]*'x^(n+1)/(n+1)) + O('x^(NM+2));
  Px = 1 - 2*B*'x + (B^2-4*C)*'x^2 + O('x^(NM+1));
  aa = lx/('x*sqrt(Px)*Fx);
  cf = vector(PMAX+1, j, polcoeff(aa, j-1));
  print("  ",NAM[k],"  a_0..a_8 = ", vector(9,j,cf[j]));
  for(i=1,#PL, my(p=PL[i], ps=psival(k,p), d);
    if(p+1>#cf, break);
    tot++;
    d = cf[p]-ps;            \\ a_{p-1} - psi(p)
    if(d%p!=0, print("    a_{p-1} != psi(p) mod p at p=",p); next);
    d = d/p;
    if((d-T[k][i])%p==0, ok1++);
    if((d+T[k][i])%p==0, ok2++);
  );
  print("    (a_{p-1}-psi(p))/p == gamma(p) mod p : ",ok1,"/",tot,"   == -gamma(p): ",ok2,"/",tot);
);
}
quit;
