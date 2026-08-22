read("apery.gp"); read("ops.gp");
default(realprecision,90);
Z3=zeta(3); G=Catalan; L32=lfun(-3,2); L33=lfun(-3,3); P2=Pi^2; P3=Pi^3; P3s=Pi^3*sqrt(3);
{ TL = ["41","46","47","49","61","84","110","111","112","133","134","135","136","142","143","228","406","~67","~88,~89","16","42","36","38","48","65","58","64","69","70"]; }
NN = if(NNarg2, NNarg2, 1200);
{
for(i=1,#OPS, my(o=OPS[i], a=o[1], hit=0);
  for(j=1,#TL, if(TL[j]==a, hit=1));
  if(!hit, next);
  my(pr=aperyPair(o[4],NN)); if(pr==0, next);
  my(A=pr[1],B=pr[2]);
  if(A[NN+1]==0, next);
  my(x=1.0*B[NN+1]/A[NN+1], xprev=1.0*B[NN-49]/A[NN-49]);
  my(ld = lindep([x,Z3,P3,P3s,P2,G,L32,L33]));
  printf("AESZ %-8s x_%d = %.30f   (drift over 50 steps: %.2e)\n", a, NN, x, abs(x-xprev));
  print("        lindep[x,z3,pi^3,pi^3sqrt3,pi^2,G,L32,L33] = ", ld~);
);
}
quit
