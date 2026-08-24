/* lattice/p2_scale/checkwin.gp -- brute-force control for the windowed
   balance search of scorel.gp `balx`.  Re-evaluates the objective
   RD^2 e_i^2 + RN^2 (q_i h22)^2 at EVERY rung of the ladder, with no window,
   and compares the argmin with the windowed answer.
   Prepend: rows_pos.gp, p2core.gp, scorel.gp.  Driver sets NLIST, KLIST, ROWF. */
{
balxfull(cf, h11, h12, h22, RD, RN) =
 my(LL=#cf, e0=h12, q0=1, e1=if(LL>=2,h12*cf[2]-h11,0), q1=if(LL>=2,cf[2],1),
    be=1, bv=RD^2*h12^2+RN^2*h22^2, vv);
 if(LL>=2, vv=RD^2*e1^2+RN^2*(q1*h22)^2; if(vv<bv, be=2; bv=vv));
 for(i=3,LL,
   my(ee=cf[i]*e1+e0, qq=cf[i]*q1+q0); e0=e1; e1=ee; q0=q1; q1=qq;
   vv = RD^2*e1^2 + RN^2*(q1*h22)^2;
   if(vv<bv, be=i; bv=vv));
 be-1;
}
PMD = ceil(7.10*vecmax(NLIST)) + 200;
default(realprecision, PMD + 30);
GG = Catalan; WD = 10^PMD; AG = floor(GG*WD);
default(realprecision, 60);
RW = rdrows(ROWF);
BADW = 0;
{
for(j=1,#NLIST,
 my(nn=NLIST[j], rw=mapget(RW,nn), XX=rw[1],YY=rw[2],VV=rw[3],UU=rw[4],
    DD=dlcm(6*nn), SS=DD^2, ZN=XX*AG-YY*WD, NN0=VV*AG-UU*WD, sZ=sign(ZN), RD=2^BB);
 for(ii=1,#KLIST,
  my(kk=KLIST[ii], TT=2^floor(kk*nn), MOD=SS*TT,
     ZHI=abs(ZN)+4*XX, NLO=NN0-4*VV, RN=(NLO*RD)\ZHI,
     hd=hermd(XX,YY,VV,UU,MOD,sZ), h11=hd[1],h12=hd[2],h22=hd[3],
     cf=contfrac(h12/h11),
     iw=balx(cf,h11,h12,h22,RD,RN)[1], ifu=balxfull(cf,h11,h12,h22,RD,RN));
  if(iw!=ifu, BADW++; print("  WINDOW MISS n=",nn," k=",kk," windowed=",iw," full=",ifu));
  print("  n=",nn," k=",kk," idx=",ifu," L=",#cf," agree=",iw==ifu)));
}
print("[checkwin] window failures: ", BADW);
\q
