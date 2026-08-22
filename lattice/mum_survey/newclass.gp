read("apery.gp"); read("ops.gp");
default(realprecision,80);
L53=lfun(5,3); L82=lfun(8,2); Z3=zeta(3); P2=Pi^2; G=Catalan;
NN=1400;
{ TL=["184","7"];
for(i=1,#OPS, my(o=OPS[i],a=o[1],h=0);
  for(j=1,#TL, if(TL[j]==a && o[3]==2, h=1));
  if(!h,next);
  my(pr=aperyPair(o[4],NN)); if(pr==0,next);
  my(A=pr[1],B=pr[2]); if(A[NN+1]==0,next);
  my(x=1.0*B[NN+1]/A[NN+1], d=abs(x-1.0*B[NN-49]/A[NN-49]));
  printf("AESZ %-5s nn=%s  x_%d = %.30f   drift=%.2e\n", a, o[2], NN, x, d);
  print("   predicted: ", if(a=="184", Str(1/4," * L(chi5,3) = ", 1.0*L53/4), Str(1/16," * L(chi8,2) = ", 1.0*L82/16)));
  print("   lindep[x, L(chi5,3), L(chi8,2), zeta(3), pi^2, G] = ", lindep([x,L53,L82,Z3,P2,G])~);
);}
quit
