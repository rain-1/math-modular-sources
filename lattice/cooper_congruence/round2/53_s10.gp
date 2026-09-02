\\ 53_s10.gp -- TASK 2a: twisted CM traces for row s10.
\\ N=10, tau_0=(3+i)/10, D0=-4, d=-4m^2, beta = 6m mod 20 (NOT 3m: 6^2+4 = 40 = 0 mod 4N),
\\ u = (eta(5t)eta(10t)/(eta(t)eta(2t)))^2, C=25, H = 1/u - 25u.
\\ H is W_5-anti-invariant, W_2-invariant  (u|W_2=u, u|W_5=1/(25u)) -- and the
\\ (W_2=+1,W_5=-1) eigenspace of L(sum of the 4 cusps) is 1-dimensional, so H is forced.
read("50_lib.gp"); read("heeg.gp");
default(realprecision, 140);
N = 10; D0 = -4; MMAX = 40;
bet = read("20_beta_s10.txt");
KAP = 5/Pi;
RR  = exp(Pi/5);
print("s10:  kappa = 1/(2 pi Im tau_0) = 5/pi = ", KAP);
print("      R = e^(pi/5) = ", RR);
print();
print("== beta-class audit: which beta mod 2N has beta^2 = d mod 4N ? ==");
{ for(m=1,6, my(d=-4*m^2, L=List());
    for(b=0,2*N-1, if((b^2-d)%(4*N)==0, listput(L,b)));
    print("  m=",m,"  d=",d,"  valid beta mod 20: ",Vec(L),"   6m mod 20 = ",(6*m)%20,
          "   3m mod 20 = ",(3*m)%20));
}
print();
{ allclass(d) = my(s=0, g=1);
  while(g^2 <= -d,
    if(d%(g^2)==0 && ((d/g^2)%4==0 || (d/g^2)%4==1) && (d/g^2)<0, s += qfbclassno(d/g^2));
    g++);
  s;
}
print("== audit A: form counts, d = -4m^2 ==");
print("  m   d      #red  sum_g h(d/g^2)   h(d)   H(-d)");
{ for(m=1,MMAX, my(d=-4*m^2);
    print("  ",m," ",d,"  ",#redforms(d),"  ",allclass(d),"  ",qfbclassno(d),"  ",qfbhclassno(4*m^2)));
}
print();
uval(t) = (eta(5*t,1)*eta(10*t,1)/(eta(t,1)*eta(2*t,1)))^2;
Hval(t) = my(w=uval(t)); 1/w - 25*w;
{ tr10(m, PB) = my(d=-4*m^2, bt=(6*m)%20, RF=redforms(d), T0=0., T1=0., nf=0, rep, al, om, ch);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,PB);
    if(rep==0, print("   *** NO REP m=",m," form=",RF[i]); next);
    nf++;
    om = omeg(rep); ch = genchar(rep,D0);
    al = (-rep[2] + I*2*m)/(2*rep[1]);
    T0 += Hval(al)/om;
    T1 += ch*Hval(al)/om);
  [T0, T1, nf, #RF];
}
print("== traces: T0 = untwisted, T1 = chi_{-4}-twisted ==");
TAU = vector(MMAX); TAU0 = vector(MMAX);
{
for(m=1,MMAX,
  my(A = tr10(m,80), t0, t1);
  t0 = A[1]; t1 = A[2];
  TAU0[m] = t0; TAU[m] = t1;
  print("m=",m,"  #found/#forms=",A[3],"/",A[4],if(A[3]==A[4],"  OK","  *** MISMATCH"));
  print("     T0 = ",t0);
  print("     T1 = ",t1,"      T1/(2i) = ",t1/(2*I));
);
}
print();
print("== integrality of T1/(2i) and T0 ==");
{ for(m=1,MMAX,
   my(a=TAU[m]/(2*I), b=TAU0[m]);
   print("  m=",m,"  T1/(2i) -> ",round(real(a)),"  err=",abs(a-round(real(a))),
         "     T0 -> ",round(real(b)),"  err=",abs(b-round(real(b)))));
}
quit;
