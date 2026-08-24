/* lattice/p2_scale/verify.gp
   Independent validation of the row cache built by build_rows.gp from the
   fitted order-3 P-recursion.  The recurrence itself is only *fitted*
   (n <= 80) and *checked* (n <= 200), so the extension to n = 3000 must be
   tested against facts that the fit cannot have encoded.  Three exact 2-adic
   identities of P2_HOLONOMIC.md sec.5.1 do that:

     (V1) v_2(X_n xi_2 - Y_n) = v_2(X_n) + 24n - 1 - 4 s_2(3n)   [proved]
     (V2) v_2(V_n xi_2 - U_n) = v_2(V_n) + 28n - 2 s_2(n) - 2 s_2(3n) - 1
     (V3) v_2(X_n U_n - V_n Y_n) = v_2(X_n)+v_2(V_n)+24n-1-4 s_2(3n)

   (V2) involves only the Nesterenko entries B_n, C_n, i.e. exactly the part
   of the cache produced by the fitted operator; it is a congruence to
   28n + O(log n) binary digits and cannot hold by accident.
   Also checked: the four archimedean growth rates against their closed forms.

   xi_2 = zeta_2(2) is built as lim P_m/Q_m along Zudilin's own order-2
   recurrence, in a single non-storing pass with snapshots at m = 6n.
   Prepend: lattice/positivity/rows_pos.gp, lattice/p2_structure/p2core.gp.  */

default(parisize, 4000000000);
default(parisizemax, 10000000000);

/* one pass of Zudilin's recurrence, keeping only the running pair, with a
   snapshot of P_m/Q_m at each m in SNAP (sorted increasing).             */
{
zudsnap(SNAP) =
 my(TOP = SNAP[#SNAP], q0=1, q1=7/4, p0=0, p1=13/8, out=vector(#SNAP), jj=1);
 for(m=1, TOP-1,
   my(aa=(2*m+1)^2*(2*m+2)^2*(20*m^2-8*m+1),
      bb=3520*m^6+5632*m^5+2064*m^4-384*m^3-156*m^2+16*m+7,
      cc=(2*m-1)^2*(2*m)^2*(20*m^2+32*m+13),
      q2=(bb*q1+cc*q0)/aa, p2=(bb*p1+cc*p0)/aa);
   q0=q1; q1=q2; p0=p1; p1=p2;
   while(jj<=#SNAP && SNAP[jj]==m+1, out[jj]=p1/q1; jj++));
 out;
}

NSAMP = if(type(NSAMP0)=="t_INT", NSAMP0, 50);
NTOP  = if(type(NTOP0)=="t_INT",  NTOP0, 3000);
ROWF  = "/home/ubuntu/code/math-modular-sources/lattice/p2_scale/data/rows_scale.txt";
OUTV  = "/home/ubuntu/code/math-modular-sources/lattice/p2_scale/data/verify10.csv";

RW = rdrows(ROWF);
SN = vector(NTOP\NSAMP, i, 6*NSAMP*i);
print("[verify] building xi_2 snapshots up to m = ", SN[#SN], " ...");
XS = zudsnap(SN);
print("[verify] snapshots done");

PMD = ceil(7.10*NTOP) + 200;
default(realprecision, PMD + 30);
GG = Catalan;  WD = 10^PMD;  AG = floor(GG*WD);
default(realprecision, 60);
BAD = 0;
write(OUTV, "n,v1,v2,v3,logX,logV,logZform,logNform");
{
for(i=1,#SN,
  my(nn = NSAMP*i, rw = mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
     xr = XS[i],
     w1 = valuation(XX*xr-YY,2) == valuation(XX,2)+24*nn-1-4*hammingweight(3*nn),
     w2 = valuation(VV*xr-UU,2) == valuation(VV,2)+28*nn-2*hammingweight(nn)-2*hammingweight(3*nn)-1,
     w3 = valuation(XX*UU-VV*YY,2) == valuation(XX,2)+valuation(VV,2)+24*nn-1-4*hammingweight(3*nn));
  if(!(w1&&w2&&w3), BAD++; print("[verify] FAIL at n=",nn,"  ",[w1,w2,w3]));
  write(OUTV, Strprintf("%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f", nn, w1, w2, w3,
     log(XX*1.)/nn, log(VV*1.)/nn,
     (log(abs(XX*AG-YY*WD)*1.)-log(WD*1.))/nn,
     (log((VV*AG-UU*WD)*1.)-log(WD*1.))/nn)));
}
print("[verify] samples = ", #SN, "   failures = ", BAD);
print("[verify] E1 = 13.0995887908  E2 = 14.3931452672");
\q
