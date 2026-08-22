default(realprecision,50);
NQ = 320;
read("lattice/weight_drop/03_setup.gp"); build();
tc = [1/27, 1/16, 1/16];
targ = [zeta(2)/7, zeta(2)/5, 0];
/* L(2,chi_-3) */
L2chi3 = sumalt(n=0, 1.0/(3*n+1)^2 - 1.0/(3*n+2)^2);
targ[3] = L2chi3/2;
{Ifun(cs,s,y0) = my(S=0.0); for(m=1,#cs, S += cs[m]*(2*Pi*m)^(-s)*incgam(s,2*Pi*m*y0)); S;}
{for(r=1,#DAT,
  my(nm=DAT[r][1],N=DAT[r][2],T=DAT[r][3],F=DAT[r][4],PHI=DAT[r][5]);
  my(y0=1/sqrt(N), q0=exp(-2*Pi*y0), M=NQ-4, cs=vector(M,m,polcoeff(PHI,m)*1.0));
  print("======== ",nm,"  N=",N,"  y0=",y0);
  /* fold test: t(i/sqrt N) = t_c */
  my(tv=subst(truncate(T),q,q0));
  print("  t(i/sqrt N) = ",tv,"   t_c = ",tc[r]*1.0,"   diff=",tv-tc[r]*1.0);
  /* Lambda(Phi,s) = I(s) - N^{2-s} I(4-s) */
  my(La(s)=Ifun(cs,s,y0)-N^(2-s)*Ifun(cs,4-s,y0));
  print("  Lambda(Phi,1)=",La(1),"  Lambda(Phi,2)=",La(2),"  Lambda(Phi,3)=",La(3));
  my(xi=4*Pi^3*La(3));
  print("  xi_pred = 4 pi^3 Lambda(Phi,3) = ",xi);
  print("  target                        = ",targ[r]*1.0);
  print("  DIFF = ",xi-targ[r]*1.0);
  /* L(Xi,s) := (2pi)^s Lambda(Xi,s)/Gamma(s),  Lambda(Xi,s)=2 pi Lambda(Phi,s+1)/s */
  my(LX(s)=(2*Pi)^s*(2*Pi*La(s+1)/s)/gamma(s));
  print("  L(Xi,1)=",LX(1),"   L(Xi,2)=",LX(2),"   L(Xi,3)=",LX(3));
  print("  L(Phi,1)=",2*Pi*La(1),"  -N/(2pi^2)*xi = ",-N/(2*Pi^2)*xi);
);}
quit;
