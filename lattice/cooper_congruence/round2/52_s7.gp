\\ 52_s7.gp -- TASK 1: independent reproduction of the s7 twisted CM trace, with
\\ class-count audits, up to m = 40.  N=7, D0=-3, d=-3m^2, beta = 5m mod 14,
\\ H = 1/u - 49u  (Fricke-anti-invariant Hauptmodul combination).
read("50_lib.gp"); read("heeg.gp");
default(realprecision, 160);
N = 7; D0 = -3; MMAX = 40;
bet = read("20_beta_s7.txt");
KAP = 7/(Pi*sqrt(3));
RR  = exp(Pi*sqrt(3)/7);
print("kappa = 1/(2 pi Im tau_0) = ", KAP);
print("R     = e^(pi sqrt3/7)    = ", RR);
print();

\\ ---- audit A: reduced-form counts vs class numbers -------------------------
\\ #redforms(d) should be sum over g with g^2|d, d/g^2 = 0,1 mod 4 of h(d/g^2)
{ allclass(d) = my(s=0, g=1);
  while(g^2 <= -d,
    if(d%(g^2)==0 && ((d/g^2)%4==0 || (d/g^2)%4==1) && (d/g^2)<0,
       s += qfbclassno(d/g^2));
    g++);
  s;
}
{ primcount(d) = my(c=0, RF=redforms(d));
  for(i=1,#RF, if(gcd(gcd(RF[i][1],RF[i][2]),RF[i][3])==1, c++)); c; }
print("== audit A: form counts, disc d = -3 m^2 ==");
print("  m   d      #red  sum_g h(d/g^2)  #prim  h(d)   sum 1/om   H(-d)");
{
for(m=1,MMAX,
  my(d=-3*m^2, RF=redforms(d), w=0.);
  for(i=1,#RF, w += 1/omeg(RF[i]));
  print("  ",m," ",d,"  ",#RF,"  ",allclass(d),"  ",primcount(d),"  ",qfbclassno(d),"  ",w,"  ",qfbhclassno(3*m^2));
);
}
print();

\\ ---- the traces -----------------------------------------------------------
uval(t) = (eta(7*t,1)/eta(t,1))^4;
Hval(t) = my(w=uval(t)); 1/w - 49*w;
{ tr7(m, PB) = my(d=-3*m^2, bt=(5*m)%14, RF=redforms(d), T=0., nf=0, rep, al, om, ch);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,PB);
    if(rep==0, print("   *** NO REP m=",m," form=",RF[i]); next);
    nf++;
    om = omeg(rep); ch = genchar(rep,D0);
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    T += ch*Hval(al)/om);
  [T, nf, #RF];
}
print("== audit B + the twisted trace tau_chi(m) = Tr/(i sqrt3) ==");
TAU = vector(MMAX);
{
for(m=1,MMAX,
  my(A = tr7(m,80), t, ti);
  t = A[1]/(I*sqrt(3));
  ti = round(real(t));
  TAU[m] = ti;
  print("m=",m,"  #found/#forms = ",A[2],"/",A[3],
        if(A[2]==A[3],"  OK","  *** MISMATCH"),
        "   tau_chi=",ti,"   |Im|=",abs(imag(t)),"   |err|=",abs(real(t)-ti));
);
}
print();
print("tau_chi(1..",MMAX,") = ", TAU);
write("52_tau_s7.txt", TAU);
print();
print("== the ratio test:  beta(m)/tau(m)  and  m beta(m)/((m-kappa) tau(m)) ==");
print("   m     beta/tau                          m*beta/((m-kappa)*tau)              resid*R^(m/2)");
{
for(m=1,MMAX,
  if(TAU[m]==0, print("   ",m,"   tau = 0  (7|m: the two beta-classes merge)"); next);
  my(r1 = bet[m]/TAU[m], r2 = m*bet[m]/((m-KAP)*TAU[m]));
  print("   ",m,"   ",r1,"   ",r2,"   ",(r2-3)*RR^(m/2));
);
}
quit;
