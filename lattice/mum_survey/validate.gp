read("apery.gp"); read("ops.gp");
default(realprecision,120);
G = Catalan; Z3=zeta(3); P2=Pi^2; Lc3 = 0;
\\ L(chi_{-3},2)
Lc3 = lfun(-3,2);
NN = 220;
targets = [["16","7/48"],["28","1/7"],["29","1/12"],["42","7/64"],["60","3/23"],["182","3/11"],["189","1/14"],["205","21/80"],["18","1/60"],["22","1/36"],["33","1/384"],["38","1/24"],["48","1/12"],["58","1/8"],["64","1/12"],["32","1/1170"]];
{
for(i=1,#OPS,
  my(o=OPS[i], a=o[1]);
  for(j=1,#targets, if(targets[j][1]==a,
    my(pr=aperyPair(o[4],NN));
    if(pr==0, print(a,": lead vanishes"); break);
    my(A=pr[1],B=pr[2], x=1.0*B[NN+1]/A[NN+1]);
    my(ld=lindep([x,1,P2,G,Lc3,Z3,Pi^4,Pi^3*sqrt(3)]));
    print("AESZ ",a,"  B/A(",NN,")=",x, "\n     lindep=",ld~, "   AvSZ says ",targets[j][2]);
  )));
}
quit
