read("apery.gp"); read("ops.gp");
default(realprecision,60);
Z3=zeta(3); G=Catalan; L32=lfun(-3,2);
{ PRED = [["41",13/108,Z3,"zeta(3)"],["46",13/54,Z3,"zeta(3)"],["47",91/864,Z3,"zeta(3)"],
   ["84",7/16,Z3,"zeta(3)"],["205",21/80,Z3,"zeta(3)"],["49",1/8,L32,"L(chi-3,2)"],
   ["61",1/648,G,"G"],["110",1/12,G,"G"],["111",1/8,G,"G"],["112",1/120,G,"G"],
   ["406",1/2,G,"G"],["133",1/8,L32,"L(chi-3,2)"],["134",1/12,L32,"L(chi-3,2)"],
   ["135",1/24,L32,"L(chi-3,2)"],["136",1/120,L32,"L(chi-3,2)"],["142",1/12,L32,"L(chi-3,2)"],
   ["143",1/120,L32,"L(chi-3,2)"],["183",-1/4,L32,"L(chi-3,2)"],["228",-15/64,L32,"L(chi-3,2)"],
   ["137",5/32,L32,"L(chi-3,2)"],["138",5/48,L32,"L(chi-3,2)"],["139",5/96,L32,"L(chi-3,2)"],
   ["140",1/96,L32,"L(chi-3,2)"]];}
NN=900;
{
for(i=1,#OPS, my(o=OPS[i], a=o[1]);
  for(j=1,#PRED, if(PRED[j][1]==a,
    my(pr=aperyPair(o[4],NN));
    if(pr==0, print(a,": skip"); break);
    my(A=pr[1],B=pr[2], x=1.0*B[NN+1]/A[NN+1], t=PRED[j][2]*PRED[j][3]);
    printf("AESZ %-6s predicted %s*%s = %.20f   measured B_%d/A_%d = %.20f   |diff| = %.3e\n",
       a, PRED[j][2], PRED[j][4], t, NN, NN, x, abs(x-t));
  )));
}
quit
