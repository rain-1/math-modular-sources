/* lattice/p2_scale/srun.gp -- one CSV line per (k,n).
   Prepend: lattice/positivity/rows_pos.gp, lattice/p2_structure/p2core.gp,
            lattice/p2_scale/scorel.gp.
   The driver must set: RW (row map), AG, WD (G = AG/WD with |G-AG/WD|<=1/WD),
   XR (the 2-adic rational representative of xi_2 = zeta_2(2)), KLIST.      */

HDR = "k,n,ok,brk,idx,parity,len,irat,pqat,pqnext,pqw5,pqw20,pqmax,logh11,logh22,logmz,logg0,incone,redst,jhi,rho,rho2,skew,loglam1,logcone,logq,v2qp,dinf,d2,q1,q2,q3,q4,q5,q6,q7,q8,q9,cY,cX,cIdx,intq";

{
anal2(nn, kk, XX, YY, VV, UU, SS, ZN, NN0, WD, XR) =
 my(TT = 2^floor(kk*nn), MOD = SS*TT, sZ = sign(ZN), RD = 2^BB);
 my(ZLO = abs(ZN)-4*XX, ZHI = abs(ZN)+4*XX, NLO = NN0-4*VV, NHI = NN0+4*VV);
 my(RNa = (NLO*RD)\ZHI, RNb = (NHI*RD)\ZLO + 1, brk = RNb-RNa);
 my(hd = hermd(XX,YY,VV,UU,MOD,sZ),
    h11=hd[1], h12=hd[2], h22=hd[3], g0=hd[4], cY=hd[6], cX=hd[7], cIdx=hd[8]);
 my(cf = contfrac(h12/h11), LL = #cf, mz = gcd(h11,h12));
 my(res = vector(2), RNv = [RNa, RNb]);
 for(s=1,2,
  my(RN = RNv[s], bl = balx(cf,h11,h12,h22,RD,RN),
     i0 = bl[1], bE = bl[2], bQ = bl[3], pE = bl[4], pQ = bl[5]);
  my(w1 = [bE, bQ*h22], w2);
  if(i0 == 0,
    /* partner: the next ladder vector */
    my(e2 = if(LL>=2, h12*cf[2]-h11, h12), q2b = if(LL>=2, cf[2], 1));
    w2 = [e2, q2b*h22],
    w2 = [pE, pQ*h22]);
  my(gr = gredx(w1, w2, RD, RN), b1 = gr[1], b2 = gr[2], st = gr[3]);
  my(q2b1 = (RD*b1[1])^2 + (RN*b1[2])^2, q2b2 = (RD*b2[1])^2 + (RN*b2[2])^2);
  my(inc = if((b1[1]>=0 && b1[2]>=0) || (b1[1]<=0 && b1[2]<=0), 1, 0));
  my(L0 = eightx(b1,b2,RD,RN), cov = RD*RN*h11*h22);
  my(jhi = if(L0==0, 8, (L0*(sqrtint(q2b1)+1))\cov + 2));
  my(sc = conex(b1,b2,RD,RN,jhi), fc = sc[1], cv = sc[2], qq2 = sc[3]);
  res[s] = [[i0, b1, b2, inc, cv, sc[4]], [st, jhi, fc, cv, qq2, q2b1, q2b2]]);
 my(same = (res[1][1] == res[2][1]), ks = res[1][1], rs = res[1][2],
    i0=ks[1], b1=ks[2], b2=ks[3], inc=ks[4],
    st=rs[1], jhi=rs[2], fc=rs[3], cv=rs[4], qq2=rs[5], q2b1=rs[6], q2b2=rs[7]);
 my(cz = sZ*cv[1], cn = cv[2],
    qn = (cz*XX+cn*VV)/MOD, pn = (cz*YY+cn*UU)/MOD,
    intq = (type(qn)=="t_INT") && (type(pn)=="t_INT"));
 my(v2qp = if(intq, valuation(qn*XR-pn, 2), 0));
 my(loglz = log(abs(ZN)*1.) - log(WD*1.) - log(MOD*1.),
    logl1 = 0.5*log(q2b1*1.) - BB*log(2.) + loglz,
    logfc = log(fc*1.) - BB*log(2.) + loglz,
    logq  = if(intq && qn, log(abs(qn)*1.), 0.));
 my(hst = vector(9), wmx5 = 0, wmx20 = 0);
 for(i=2,LL, my(a=cf[i]); if(a<=8, hst[a]++, hst[9]++));
 for(i=max(2,i0+2-5), min(LL,i0+2+5), if(cf[i]>wmx5, wmx5=cf[i]));
 for(i=max(2,i0+2-20), min(LL,i0+2+20), if(cf[i]>wmx20, wmx20=cf[i]));
 write(OUTF, Strprintf("%.4f,%d,%d,%d,%d,%d,%d,%.6f,%d,%d,%d,%d,%d,%.5f,%.5f,%.5f,%.5f,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
   kk, nn, if(same && cY && cX && cIdx && intq && (sZ==(-1)^nn) && brk<=2, 1, 0),
   brk, i0, i0%2, LL, if(LL>1, i0/(LL-1)*1., 0.),
   if(i0+1<=LL, cf[i0+1], 0), if(i0+2<=LL, cf[i0+2], 0), wmx5, wmx20, vecmax(cf),
   log(h11*1.)/nn, log(h22*1.)/nn, log(mz*1.)/nn, log(g0*1.)/nn,
   inc, st, jhi,
   fc/sqrt(q2b1*1.), sqrt(qq2*1./q2b1), sqrt(q2b2*1./q2b1),
   logl1/nn, logfc/nn, logq/nn, v2qp,
   if(logq, -logfc/logq, 0.), if(logq, v2qp*log(2.)/logq, 0.),
   hst[1],hst[2],hst[3],hst[4],hst[5],hst[6],hst[7],hst[8],hst[9],
   cY, cX, cIdx, if(intq,1,0)));
}
