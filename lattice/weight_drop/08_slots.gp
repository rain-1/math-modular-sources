default(parisizemax,4000000000);
default(realprecision,60);
NQ = 420;
read("lattice/weight_drop/03_setup.gp"); build();
{Ifun(cs,s,y0)=my(S=0.0); for(m=1,#cs, S += cs[m]*(2*Pi*m)^(-s)*incgam(s,2*Pi*m*y0)); S;}
Lchi7_1 = lfun(-7,1); Lchi3_3 = lfun(-3,3); Lchi3_2=lfun(-3,2);
{for(r=1,#DAT, my(nm=DAT[r][1],N=DAT[r][2],PHI=DAT[r][5]);
  my(y0=1/sqrt(N), M=NQ-4, cs=vector(M,m,polcoeff(PHI,m)*1.0));
  my(La(s)=Ifun(cs,s,y0)-N^(2-s)*Ifun(cs,4-s,y0));
  my(LX(s)=(2*Pi)^s*(2*Pi*La(s+1)/s)/gamma(s));
  print("=== ",nm);
  print("  Lambda(Phi,1)=L(Phi,1)/2pi ; L(Phi,1) = ",2*Pi*La(1));
  print("  L(Xi,2)=xi = ",LX(2));
  print("  L(Xi,3)   = ",LX(3));
  print("  L(Xi,4)   = ",LX(4));
);}
print();
print("s7  : L(Xi,3)/(zeta(3)*L(1,chi-7)) = ", 0);
quit;
