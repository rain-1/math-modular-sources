\\ 59_s18ident.gp -- TASK (B): identification of beta_{s18} as a twisted CM trace of
\\ fhat = D f + f/(2 pi y),  f = 1/(x F)  on Gamma_0(18).
\\ Two fixes to the earlier attempt:
\\  (1) the admissible Heegner set is {Q of disc -36m^2 : 3 does not divide cont(Q)} --
\\      forms with 3 | cont have NO Gamma_0(18)-Heegner representative at all, because
\\      Q = 3Q' needs 6 | Q'(p,r) and 3 | p^2+r^2 forces 3|p, 3|r, contra gcd(p,r)=1. [proved]
\\  (2) fhat is computed by BOTH routes and the q-series route is used whenever the
\\      closed form is unsafe (|F| small OR the two disagree).
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N = 18; MMAX = 22;
initfser(3,900);
bet = read("20_beta_s18.txt");
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
\\ robust fhat: q-series route, cross-checked against the closed form when the latter is safe
{ fh(t) = my(Fv=Fmod(3,t), vq, vc);
  vq = fhatQ(t, 18, FSER[3], #Vec(FSER[3][1])-3);
  if(abs(Fv) > 1e-12,
     vc = fhatC(3,t)[2];
     if(abs(vc-vq) > 1e-20*(1+abs(vq)), print("   !! fhat mismatch at ",t,": ",abs(vc-vq))));
  vq;
}
{ tr18(m, bt, PB) = my(d=-36*m^2, RF=redforms(d), T=[0.,0.,0.], nf=0, nadm=0, rep, al, om, h);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    nadm++;
    rep = heegrep(RF[i],N,bt,PB);
    if(rep==0, print("   *** NO REP m=",m," ",RF[i]); next);
    nf++;
    om = omeg(rep);
    al = (-rep[2] + I*6*m)/(2*rep[1]);
    h = fh(al)/om;
    T[1] += h;
    T[2] += genchar(rep,-4)*h;
    T[3] += genchar(rep,-3)*h);
  [T, nf, nadm, #RF];
}
print("== s18: beta(m) vs the three traces, beta-class 18m mod 36 ==");
print("   lambda candidates: -nu^2/g'(u_0) = -16/(8 sqrt3) = ", -16/(8*sqrt(3)));
{
for(m=1,MMAX,
  my(A=tr18(m,(18*m)%36,80), T=A[1]);
  print("m=",m,"  #found/#adm/#all=",A[2],"/",A[3],"/",A[4],if(A[2]==A[3],"","  ***MISSING***"),"   beta=",bet[m]);
  print("     triv  : ",T[1],"    beta/T = ",if(abs(T[1])>1e-30, bet[m]/T[1], "-"));
  print("     chi_-4: ",T[2],"    beta/T = ",if(abs(T[2])>1e-30, bet[m]/T[2], "-"));
  print("     chi_-3: ",T[3],"    beta/T = ",if(abs(T[3])>1e-30, bet[m]/T[3], "-"));
);
}
print();
print("== all admissible beta mod 36, chi_-3 and trivial, m=1..8 ==");
{
for(m=1,8,
  my(d=-36*m^2);
  print("--- m=",m,"  beta_row=",bet[m]);
  for(b=0,35,
    if((b^2-d)%(4*N)!=0, next);
    my(A=tr18(m,b,80), T=A[1]);
    print("    beta=",b,"  found ",A[2],"/",A[3],
          "   triv=",T[1],"   chi_-3=",T[3],"   chi_-4=",T[2]));
);
}
quit;
