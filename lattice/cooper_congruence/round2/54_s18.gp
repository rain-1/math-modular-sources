\\ 54_s18.gp (v2) -- TASK 2b: CM traces for row s18.
\\ N=18, tau_0=(3+i)/6 = fixed point of W_9, form [18,-18,5], disc -36.
\\ d = -36 m^2, beta = 18m mod 36  (forced by the phase: q_0^{-m}=(-1)^m R^m needs
\\ exp(2 pi i beta/36) = (-1)^m).  H = 1/u - u.
\\ NOTE: forms with 3 | content have NO Gamma_0(18)-Heegner representative
\\ (Q = 3Q' needs 6 | Q'(p,r), impossible since 3 | p^2+r^2 forces 3|p,3|r).
\\ So the admissible set is {Q : 3 does not divide content(Q)}; we check that ALL of
\\ those are found.
read("50_lib.gp"); read("heeg.gp");
default(realprecision, 140);
N = 18; MMAX = 30;
bet = read("20_beta_s18.txt");
KAP = 3/Pi;  RR = exp(Pi/3);
print("s18:  kappa = 3/pi = ", KAP, "     R = e^(pi/3) = ", RR);
print();
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
uval(t) = (eta(2*t,1)*eta(3*t,1)^2*eta(18*t,1)/(eta(t,1)*eta(6*t,1)^2*eta(9*t,1)))^6;
Hval(t) = my(w=uval(t)); 1/w - w;
{ tr18(m, PB) = my(d=-36*m^2, bt=(18*m)%36, RF=redforms(d), T=[0.,0.,0.,0.], nf=0, nadm=0, rep, al, om, h);
  for(i=1,#RF,
    if(cont(RF[i])%3==0, next);
    nadm++;
    rep = heegrep(RF[i],N,bt,PB);
    if(rep==0, print("   *** NO REP m=",m," form=",RF[i]); next);
    nf++;
    om = omeg(rep);
    al = (-rep[2] + I*6*m)/(2*rep[1]);
    h = Hval(al)/om;
    T[1] += h;
    T[2] += genchar(rep,-4)*h;
    T[3] += genchar(rep,-3)*h;
    if(cont(RF[i])==1, T[4] += genchar(rep,-3)*h));
  [T, nf, nadm, #RF];
}
print("== traces (admissible forms = content coprime to 3) ==");
V=vector(MMAX); VP=vector(MMAX); V1=vector(MMAX); V2=vector(MMAX);
{
for(m=1,MMAX,
  my(A=tr18(m,80), T=A[1]);
  V1[m]=T[1]; V2[m]=T[2]; V[m]=T[3]; VP[m]=T[4];
  print("m=",m,"  #found/#admissible/#all = ",A[2],"/",A[3],"/",A[4],
        if(A[2]==A[3],"  OK","  *** MISSING REPS"));
  print("    triv        : ",T[1]);
  print("    chi_-4      : ",T[2]);
  print("    chi_-3      : ",T[3],"     /sqrt(3) = ",T[3]/sqrt(3));
  print("    chi_-3 prim : ",T[4],"     /sqrt(3) = ",T[4]/sqrt(3));
);
}
print();
print("== integrality of T_{chi_-3}/sqrt(3) ==");
TT = vector(MMAX);
{ for(m=1,MMAX,
   my(a = V[m]/sqrt(3), r);
   r = round(real(a));
   if(abs(a-r) < 1e-60, TT[m]=r);
   print("  m=",m,"  T/sqrt3 = ",a,"     -> integer? ",if(abs(a-r)<1e-60,concat("YES ",Str(r)),"NO")));
}
print();
print("tau_chi(1..",MMAX,") = ",TT,"   (0 = not an integer / anomalous)");
write("54_tau_s18.txt", TT);
print();
print("== ratio test (odd m):  beta/tau  and  m beta/((m-kappa) tau) ==");
{ for(m=1,MMAX,
   if(TT[m]==0, print("   ",m,"   anomalous (beta = 0 mod 36, the +-beta classes merge)"); next);
   my(r1=bet[m]/TT[m], r2=m*bet[m]/((m-KAP)*TT[m]));
   print("   m=",m,"   beta/tau=",r1,"   r2=",r2));
}
print();
print("== residual (r2 - L)*R^(m/2) for L = 3, 4, 6, 8, 12 ==");
{ for(m=1,MMAX,
   if(TT[m]==0, next);
   my(r2=m*bet[m]/((m-KAP)*TT[m]));
   print("   m=",m,"  ",vector(5,i,my(L=[3,4,6,8,12][i]); (r2-L)*RR^(m/2))));
}
quit;
