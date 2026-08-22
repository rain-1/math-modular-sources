/* ================================================================
   survey.gp -- Apery limits of order-4 Calabi-Yau operators, at the
   archimedean place and at every small prime.

   For each AESZ operator (data from CYCluster, ops.gp):
     A_n, B_n   the Apery pair (Almkvist-van Straten-Zudilin)
     lam1,lam2  growth of A_n and of the linear form B_n - xi A_n
     k          minimal exponent with d_n^k B_n in Z
     xi         archimedean Apery limit, identified against
                zeta(3), pi^2, G, L(chi-3,2), L(chi-3,3), pi^4, Q
     sigma_p    p-adic slope of B_n/A_n
     xi_p       p-adic Apery limit, identified against Kubota-Leopoldt
                values L_p(s, chi omega^m).
   Output: one "|"-separated line per operator.
   ================================================================ */
read("apery.gp"); read("ops.gp"); read("lpgen.gp");
default(realprecision, 260);
NN  = if(NNarg, NNarg, 200);
PRIMES = [2,3,5,7,11,13];
{CONS = [["Q",1],["zeta(3)",zeta(3)],["pi^2",Pi^2],["G",Catalan],
        ["L(chi-3,2)",lfun(-3,2)],["L(chi-3,3)",lfun(-3,3)],["pi^4",Pi^4]];}

idArch(x) =
{ for(j=1,#CONS,
    my(r = x/CONS[j][2], q = bestappr(r, 10^8));
    if(q!=0 && abs(r-q) < abs(r)*10^(-40), return([CONS[j][1], q])));
  ["?", 0];
}
padicTargets(p, PR) =
 {[["zeta_p(3)",        LpG(p,triv, -2,3,PR)],
  ["L_p(2,om^-1)",     LpG(p,triv, -1,2,PR)],
  ["L_p(2,chi-4om^-1)",LpG(p,chim4,-1,2,PR)],
  ["L_p(2,chi-3om^-1)",LpG(p,chim3,-1,2,PR)],
  ["L_p(3,chi-3om^-2)",LpG(p,chim3,-2,3,PR)],
  ["L_p(3,chi-4om^-2)",LpG(p,chim4,-2,3,PR)]];}
{
print("#aesz|nn|degz|k|lam1|lam2|archconst|archrat|slopes[p,sigma,vp(A_N)]|padic[p,digits,const,rat]");
for(i=1,#OPS,
  my(o=OPS[i], a=o[1], nn=o[2], dz=o[3], pols=o[4]);
  my(pr = aperyPair(pols, NN));
  if(pr==0, print(a,"|",nn,"|",dz,"|SKIP"); next);
  my(A=pr[1], B=pr[2]);
  if(A[NN+1]==0, print(a,"|",nn,"|",dz,"|SKIP"); next);
  my(x = 1.0*B[NN+1]/A[NN+1], idr = idArch(x));
  my(lam1 = exp(log(abs(1.0*A[NN+1]))/NN));
  my(m = NN\2, lin = abs(1.0*B[m+1] - x*A[m+1]), lam2 = if(lin==0, 0., exp(log(lin)/m)));
  my(kk = denomExp(B, min(NN,120)));
  my(sl=[], pd=[]);
  for(t=1,#PRIMES,
    my(p=PRIMES[t], n1=NN\2, n2=NN);
    my(v1=slopeVal(A,B,p,n1), v2=slopeVal(A,B,p,n2));
    my(sig = (v2-v1)*1.0/(n2-n1), ka = -valuation(A[NN+1],p)*1.0/NN);
    sl = concat(sl, [[p, round(sig*100)/100., round(ka*100)/100.]]);
    if(sig > 0.2,
      my(PRp = max(10, floor(sig*NN*0.85)));
      my(xi = B[NN+1]/A[NN+1] + O(p^PRp), TG = padicTargets(p,PRp), hits=[]);
      for(u=1,#TG,
        if(TG[u][2]==0, next);
        my(rel = lindep([xi, TG[u][2]]));
        if(#rel==2 && rel[1]!=0 && abs(rel[1])<10^7 && abs(rel[2])<10^7,
           if(valuation(xi*rel[1]+TG[u][2]*rel[2],p) > PRp-10,
              hits = concat(hits,[[TG[u][1], -rel[2]/rel[1]]]))));
      if(#hits==0 && valuation(xi,p) > PRp-10, hits=[["ZERO",0]]);
      pd = concat(pd, [[p, PRp, hits]])));
  print(a,"|",nn,"|",dz,"|",kk,"|",lam1,"|",lam2,"|",idr[1],"|",idr[2],"|",sl,"|",pd);
);
}
quit
