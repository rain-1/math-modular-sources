/* lattice/p2_scale/run_scan.gp -- driver.
   The caller must define, before this file:  NL, NH, KLIST, OUTF, ROWF,
   and optionally NSTEP (default 1; = NCH when the split is interleaved).
   cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
       lattice/p2_scale/scorel.gp lattice/p2_scale/srun.gp \
       PARAMS lattice/p2_scale/run_scan.gp > run.gp && gp -q run.gp          */

PMD = ceil(7.10*NH) + 200;   /* decimal digits of G kept exactly: must beat
                                 both horizons, 6.492 n (Nesterenko) plus the
                                 0.562 n growth of r = lambda_N/lambda_Z      */
default(realprecision, PMD + 30);
GG  = Catalan;
WD  = 10^PMD;
AG  = floor(GG*WD);
default(realprecision, 60);
MTOP = NH + 200;                    /* 2-adic precision of xi_2: > 8*MTOP-1 */
ZZ  = zud(MTOP);
XR  = ZZ[2][MTOP+1]/ZZ[1][MTOP+1];
ZZ  = 0;
RW  = rdrows(ROWF);
write(OUTF, HDR);
{
forstep(nn = NL, NH, if(type(NSTEP)=="t_INT", NSTEP, 1),
  my(rw = mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
     DD = dlcm(6*nn), SS = DD^2, ZN = XX*AG - YY*WD, NN0 = VV*AG - UU*WD);
  for(ii=1,#KLIST, anal2(nn, KLIST[ii], XX,YY,VV,UU, SS, ZN, NN0, WD, XR)));
}
\q
