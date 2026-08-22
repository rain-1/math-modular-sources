default(realprecision, 340);
/* synthetic: does the pair-lindep in b8's style detect a planted relation? */
x = 3/7*zeta(3) - 11/13*lfun(-52,4);
r = lindep([x, zeta(3), lfun(-52,4)], 300);
print("planted 2-term: ", r~);
y = 1/5 + 2/9*Pi^4 + 7/3*lfun(17,3);
print("planted 3-term: ", lindep([y,1,Pi^4,lfun(17,3)],290)~);
/* pipeline validation on AESZ 16 (known 7/48 zeta(3)) */
default(parisize,4000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
{ my(nms=["16","42","28","29","182","185","58"]);
  for(i=1,#OPS, my(o=OPS[i], h=0);
    for(j=1,#nms, if(nms[j]==o[1], h=1));
    if(!h, next);
    my(pr=aperyPair(o[4],600)); if(pr==0,next);
    my(A=pr[1],B=pr[2]); if(A[601]==0,next);
    my(x=1.0*B[601]/A[601], d=x-1.0*B[600]/A[600]);
    printf("AESZ %-4s nn=%-8s xi=%.25f  drift=%.2e  ", o[1], o[2], x, d);
    my(hit=0);
    foreach([["zeta(3)",zeta(3)],["Pi^2",Pi^2],["G",Catalan],["L(-3,2)",lfun(-3,2)],["L(-3,3)",lfun(-3,3)]], C,
      my(q=bestappr(x/C[2],10^8)); if(q!=0 && abs(x/C[2]-q)<1e-30, print("= ",q," * ",C[1]); hit=1));
    if(!hit, print("?"))); }
quit
