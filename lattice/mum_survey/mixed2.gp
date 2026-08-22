read("apery.gp"); read("ops.gp");
default(realprecision,120);
Z3=zeta(3); G=Catalan; L32=lfun(-3,2); P2=Pi^2; P3=Pi^3; P3s=Pi^3*sqrt(3);
/* [aesz, predicted even part, even constant, family]  fam: 3 = zeta(3), 2G, 2L */
{ TL = [["41",13/108,3],["46",13/54,3],["47",91/864,3],["84",7/16,3],["205",21/80,3],
        ["49",1/8,21],["133",1/8,21],["134",1/12,21],["135",1/24,21],["136",1/120,21],
        ["142",1/12,21],["143",1/120,21],["183",-1/4,21],["228",-15/64,21],["137",5/32,21],
        ["61",1/648,22],["110",1/12,22],["111",1/8,22],["112",1/120,22],["406",1/2,22],
        ["~67",5/216,22],["~88,~89",1/24,22]]; }
NN=1400;
{
for(i=1,#OPS, my(o=OPS[i], a=o[1], idx=0);
  for(j=1,#TL, if(TL[j][1]==a, idx=j));
  if(!idx, next);
  my(pr=aperyPair(o[4],NN)); if(pr==0, next);
  my(A=pr[1],B=pr[2]); if(A[NN+1]==0, next);
  my(x=1.0*B[NN+1]/A[NN+1], drift=abs(x-1.0*B[NN-49]/A[NN-49]));
  my(fam=TL[idx][3], ev=TL[idx][2]);
  my(base = if(fam==3, [Z3,P3s,P3], if(fam==21, [L32,P2,P3], [G,P2,P3])));
  my(rem = x - ev*base[1]);
  my(dig = if(drift==0, 100, floor(-log(drift)/log(10))));
  printf("AESZ %-9s x=%.16f  accurate~1e-%d   even part %s*const\n", a, x, dig, ev);
  if(dig>=6,
    default(realprecision, max(20, dig-1));
    my(r2 = rem*1.0);
    my(l1=lindep([r2,base[2]]), l2=lindep([r2,base[3]]), l3=lindep([r2,base[2],base[3]]));
    printf("     remainder %.14f :  vs c2: %s   vs c3: %s   vs both: %s\n", r2, l1~, l2~, l3~);
    default(realprecision,120));
);
}
quit
