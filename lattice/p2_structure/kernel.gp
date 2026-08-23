/* lattice/p2_structure/kernel.gp   --- Task 2: the rationality-kernel vector.
   Prepend rows_pos.gp, p2core.gp.

   If G = a/b exactly then the direction
        c^ker = (a V_n - b U_n, -(a X_n - b Y_n))
   annihilates the form: c^ker.(X G - Y, V G - U) = 0 identically, and its
   two scaled coordinates are equal in absolute value, namely
        |c^ker_Z| * lz = |c^ker_N| * ln = b * MOD * lz * ln,
   so its scaled length is  sqrt(2) * b * MOD * lz * ln = sqrt(2) * b *
   (MOD/idx) * covol.  Consequences (all exact):
     -- the kernel vector beats the first minimum iff b < b_crit :=
        lambda_1 / (sqrt(2) (MOD/idx) covol),  and (1/n) log b_crit -> |F|;
     -- since lambda_1 is *measured* with the true G, the same inequality
        run backwards is a (weak) rigorous statement: G is not a rational
        of denominator <= b_crit(n).
   kerchk() verifies the identity for a rational surrogate; bcrit() tabulates
   the critical denominator in the true-G metric.                            */

/* exact check of the kernel identity for the surrogate GS = a/b */
{
kerchk(n, klist, GS, ee, RW) =
 my(rw=mapget(RW,n), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*n,i,i)), SS=DD^2,
    aa=numerator(GS), bb=denominator(GS));
 for(ii=1,#klist,
  my(kk=klist[ii], TT=2^floor(kk*n), MOD=SS*TT);
  my(B0=kfull(XX,YY,VV,UU,MOD), idx=abs(matdet(B0)));
  my(LZ=(XX*GS-YY)/MOD, LN=(VV*GS-UU)/MOD, lz=abs(LZ), ln=abs(LN));
  my(ck=[aa*VV-bb*UU, -(aa*XX-bb*YY)], g=gcd(ck[1],ck[2]), c0=[ck[1]/g,ck[2]/g]);
  my(A1=c0[1]*XX+c0[2]*VV, A2=c0[1]*YY+c0[2]*UU,
     tt=MOD/gcd(MOD,gcd(A1,A2)), cm=[tt*c0[1], tt*c0[2]]);
  my(zerof = (cm[1]*LZ+cm[2]*LN)==0,          /* form vanishes exactly */
     equal = (abs(ck[1])*lz)==(abs(ck[2])*ln), /* the two coords agree */
     pred  = bb*MOD*lz*ln,                     /* = |c^ker_Z| lz          */
     inlat = (denominator((cm[1]*XX+cm[2]*VV)/MOD)==1)
           && (denominator((cm[1]*YY+cm[2]*UU)/MOD)==1),
     kerlen= sqrt(2.)*(1.*tt/g)*pred,
     cov   = idx*lz*ln*1.);
  printf("%d,%.4f,%d,%d,%d,%d,%.6f,%.6f,%d,%d,%.6f,%.6f\n",
    ee, kk, n, zerof, equal, inlat,
    log(abs(ck[1])*lz*1.)/n, log(1.*pred)/n, tt, g,
    log(kerlen)/n, (log(1.*bb)+log(cov))/n + log(2.)/(2*n)));
}

/* critical denominator in the true-G metric */
{
bcrit(n, klist, GG, RW) =
 my(rw=mapget(RW,n), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*n,i,i)), SS=DD^2);
 for(ii=1,#klist,
  my(kk=klist[ii], TT=2^floor(kk*n), MOD=SS*TT);
  my(B0=kfull(XX,YY,VV,UU,MOD), idx=abs(matdet(B0)));
  my(LZ=(XX*GG-YY)/MOD, LN=(VV*GG-UU)/MOD,
     sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN));
  my(BB=[sZ*B0[1,1], sZ*B0[1,2]; sN*B0[2,1], sN*B0[2,2]]);
  my(gr=gred(BB,lz,ln), l1=gr[2], l2=gr[3], cov=idx*lz*ln);
  my(bc = l1/(sqrt(2.)*(MOD/idx)*cov));
  printf("%.4f,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n",
    kk, n, log(bc)/n, log(cov)/(2*n), log(l1)/n, log(l2)/n, log(1.*MOD/idx)/n));
}
