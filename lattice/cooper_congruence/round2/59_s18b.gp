\\ 59_s18b.gp -- s18 identification, with the two fixes:
\\   (1) admissible forms = {3 does not divide cont(Q)}  [the rest have NO Heegner rep, proved]
\\   (2) SMALLEST-A Heegner representative -> maximal Im -> no numerical blow-up.
\\ (The blow-ups in the earlier run came from heegrep returning the FIRST rep in its scan,
\\  e.g. A = 5634 instead of A = 18 for [1,0,9]: Im = 3m/A, so the q-series diverged and the
\\  closed form lost all precision.)
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N=18; initfser(3,900);
bet = read("20_beta_s18.txt");
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
\\ robust fhat: closed form unless F is small, then q-series; report disagreement
{ fh(t) = my(Fv=Fmod(3,t), vq, vc);
  if(abs(Fv) < 1e-15, return(fhatQ(t,18,FSER[3],#Vec(FSER[3][1])-3)));
  vc = fhatC(3,t)[2];
  vq = fhatQ(t,18,FSER[3],#Vec(FSER[3][1])-3);
  if(abs(vc-vq) > 1e-25*(1+abs(vq)), print("      [warn] fhat routes differ at Im=",imag(t),": ",abs(vc-vq)));
  vc;
}
LAM = -16/(8*sqrt(3));
print("s18:  lambda candidate -nu^2/g'(u_0) = ", LAM, "   (|lambda| = 2/sqrt3 = ",2/sqrt(3),")");
print();
{
for(m=1,20,
  my(d=-36*m^2, bt=(18*m)%36, RF=redforms(d), T=[0.,0.,0.], nadm=0, nf=0, Amax=0, Amin=10^9);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    nadm++;
    my(rep=heegmin(RF[i],N,bt,60), om, al, h);
    if(rep==0, print("   NOREP ",RF[i]); next);
    nf++;
    Amax = max(Amax, rep[1]); Amin = min(Amin, rep[1]);
    om = omeg(rep);
    al = (-rep[2] + I*6*m)/(2*rep[1]);
    h = fh(al)/om;
    T[1] += h; T[2] += genchar(rep,-4)*h; T[3] += genchar(rep,-3)*h);
  print("m=",m,"  adm/all=",nadm,"/",#RF,"  found=",nf,"  A in [",Amin,",",Amax,"]   beta=",bet[m]);
  print("    triv   : ",T[1],"   beta/T=",if(abs(T[1])>1e-25,bet[m]/T[1],"-"));
  print("    chi_-4 : ",T[2],"   beta/T=",if(abs(T[2])>1e-25,bet[m]/T[2],"-"));
  print("    chi_-3 : ",T[3],"   beta/T=",if(abs(T[3])>1e-25,bet[m]/T[3],"-"));
);
}
quit;
