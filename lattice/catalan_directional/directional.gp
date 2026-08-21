/* Directional-mass audit (Sol's proposal).  Run: cat rows_common.gp directional.gp > run.gp; gp -q run.gp */
\p 60
{
report(nlist) = my();
 printf("n  |  logPhi_rig/n  logPhi_opt/n  logPhiNET_rig/n  logPhiNET_opt/n  | symN_opt/n symNET_N/n | logS/n\n");
 for(i=1,#nlist,
  my(n=nlist[i], D=lcm(vector(6*n,i2,i2)), S=D^2,
     zr=zudrow(n), nr=nestrow(n), X=zr[1],Y=zr[2],V=nr[1],U=nr[2],
     LR=0.0, LO=0.0, LNR=0.0, LNO=0.0, MO=0.0, MNO=0.0, wl=0.0, forcedtot=0.0);
  forprime(l=3, 6*n,
   my(vX=valuation(X,l), vY=valuation(Y,l), vV=valuation(V,l), vU=valuation(U,l),
      s=2*valuation(D,l), L=log(l));
   my(mZ=min(vX,vY), mN=min(vV,vU));
   my(eO=max(0,mN-mZ), eR=if(vX!=vY, eO, 0));
   my(fZ=max(0, min(vU,s)-vY), fN=max(0, min(vY,s)-vU));
   LR += eR*L; LO += eO*L;
   LNR += max(0,eR-fZ)*L; LNO += max(0,eO-fZ)*L;
   my(gO=max(0,mZ-mN), gR=if(vV!=vU, gO, 0));
   MO += gR*L; MNO += max(0,gR-fN)*L;
   forcedtot += fZ*L;
   if(l>4*n && l<6*n, wl += eO*L);
  );
  printf("%3d | %12.5f %12.5f %14.5f %14.5f | %10.5f %10.5f | %7.4f  (forced_cZ/n=%.4f, window(4n,6n) opt=%.4f)\n",
     n, LR/n, LO/n, LNR/n, LNO/n, MO/n, MNO/n, log(S)/n, forcedtot/n, wl/n);
 );
}
report([5,6,8,10,12,15,18,20,24,28,32,36,40,46,52,60]);
/* exponent bookkeeping */
{
 my(ph=(1+sqrt(5))/2, Tp=(3303+437*sqrt(57))/144, Tm=(3303-437*sqrt(57))/144);
 my(A1=12*log(2)+12+15*log(ph), E1=12*log(2)+12-15*log(ph),
    A2=14*log(2)+12+2*log(Tp), E2=14*log(2)+12+2*log(Tm));
 my(ks=(E1+E2-12)/log(2));
 printf("\nA1=%.10f E1=%.10f A2=%.10f E2=%.10f  k*=%.10f\n",A1,E1,A2,E2,ks);
 printf(" k     sigma     F(k)      H(k)      beta_Z=x    beta_N=kappa-x   delta=H/(H-F)... \n");
 for(j=0,4, my(k=[ks,22.4,22.8,23.0,23.9][j+1], sig=12+k*log(2), x=(sig+E2-E1)/2,
     F=x+E1-sig, H=sig-x+A2-sig);
   H=F+A2-E2;
   printf("%6.3f %9.5f %9.5f %9.5f %11.5f %11.5f\n", k, sig, F, H, x, sig-x));
}
\q
