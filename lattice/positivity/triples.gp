/* lattice/positivity/triples.gp --- Task 3 (c): does a THIRD positive row help?
   Same setup as pairs.gp but with three moments (m,j1),(m,j2),(m,j3):
   MOD = d * 2^{v}, v = v_2(gcd of the three 2x2 mixed minors);
   K = {p in Z^3 : MOD | p.a, MOD | p.b}  (rank 3, index <= MOD^2);
   cone minimum of p1*M1+p2*M2+p3*M3 over K cap Z^3_{>=0} \ {0}.
   Prepend rows_pos.gp, cone80.gp, pairs.gp.                                 */

{
klat3(a1,b1,a2,b2,a3,b3,MOD) =
 my(K=matkerint([a1,a2,a3,MOD,0,0,0,0; b1,b2,b3,0,MOD,0,0,0]));
 mathnf(matrix(3,#K[1,],i,j,K[i,j]));
}

{
tripstat(AA,BB,MM, jj, RR) =
 localprec(700);
 my(A1=AA[jj[1]+1], A2=AA[jj[2]+1], A3=AA[jj[3]+1],
    B1=BB[jj[1]+1], B2=BB[jj[2]+1], B3=BB[jj[3]+1],
    M1=MM[jj[1]+1], M2=MM[jj[2]+1], M3=MM[jj[3]+1]);
 my(dd=lcm(lcm(lcm(denominator(A1),denominator(B1)),
               lcm(denominator(A2),denominator(B2))),
           lcm(denominator(A3),denominator(B3))));
 my(a1=dd*A1,b1=dd*B1,a2=dd*A2,b2=dd*B2,a3=dd*A3,b3=dd*B3);
 my(gg=gcd(gcd(a1*b2-a2*b1, a1*b3-a3*b1), a2*b3-a3*b2));
 my(v2h=if(gg==0,0,valuation(gg,2)), MOD=dd*2^v2h);
 my(B0=klat3(a1,b1,a2,b2,a3,b3,MOD), idx=abs(matdet(B0)));
 my(l1=dd*M1/MOD, l2=dd*M2/MOD, l3=dd*M3/MOD);
 /* LLL in the metric diag(l1,l2,l3) */
 my(sc=10^(150+ceil(log(1/vecmin([l1,l2,l3]))/log(10))));
 my(GR=[round(B0[1,1]*l1*sc),round(B0[1,2]*l1*sc),round(B0[1,3]*l1*sc);
        round(B0[2,1]*l2*sc),round(B0[2,2]*l2*sc),round(B0[2,3]*l2*sc);
        round(B0[3,1]*l3*sc),round(B0[3,2]*l3*sc),round(B0[3,3]*l3*sc)]);
 my(Br=B0*qflll(GR));
 my(best=0, bv=0);
 for(i=-RR,RR, for(j=-RR,RR, for(k=-RR,RR,
   my(c1=i*Br[1,1]+j*Br[1,2]+k*Br[1,3],
      c2=i*Br[2,1]+j*Br[2,2]+k*Br[2,3],
      c3=i*Br[3,1]+j*Br[3,2]+k*Br[3,3]);
   if(c1>=0 && c2>=0 && c3>=0 && (c1||c2||c3),
     my(v=c1*l1+c2*l2+c3*l3);
     if(best==0 || v<best, best=v; bv=[c1,c2,c3])))));
 if(best==0, return(0));
 my(q=(bv[1]*a1+bv[2]*a2+bv[3]*a3)/MOD, p=(bv[1]*b1+bv[2]*b2+bv[3]*b3)/MOD);
 my(ok=(type(q)=="t_INT")&&(type(p)=="t_INT")&&(q!=0)
      &&abs(abs(q*Catalan-p)-best)<best*10.^(-precision(1.*l1)+50));
 [best, if(q==0,0,log(abs(q))), if(ok,1,0), v2h, log(idx)];
}

{
runtrip(m, RR) =
 my(NN=m+1, AA=vector(NN), BB=vector(NN), MM=vector(NN));
 for(j=0,m, my(r=mom(m,j)); AA[j+1]=r[1]; BB[j+1]=r[2]; MM[j+1]=r[1]*Catalan-r[2]);
 /* best PAIR at this m, for reference */
 my(bp=0, bpj=0, bpq=0);
 for(j1=0,m-1, for(j2=j1+1,m,
   my(s=pairstat(AA[j1+1],BB[j1+1],AA[j2+1],BB[j2+1],MM[j1+1],MM[j2+1],1,""));
   if(s[5]==1 && (bp==0 || s[1]<bp), bp=s[1]; bpj=[j1,j2]; bpq=s[4])));
 printf("m=%d  best PAIR   (%d,%d): rate/m=%.5f  logq/m=%.4f  delta=%.4f\n",
   m, bpj[1], bpj[2], log(bp)/m, bpq/m, -log(bp)/bpq);
 /* triples: a window around the best pair */
 my(bt=0, btj=0, btq=0, lo=max(0,bpj[1]-4), hi=min(m,bpj[2]+4));
 for(j1=lo,hi-2, for(j2=j1+1,hi-1, for(j3=j2+1,hi,
   my(s=tripstat(AA,BB,MM,[j1,j2,j3],RR));
   if(type(s)=="t_VEC" && s[3]==1 && (bt==0 || s[1]<bt), bt=s[1]; btj=[j1,j2,j3]; btq=s[2]))));
 if(bt!=0,
   printf("m=%d  best TRIPLE (%d,%d,%d): rate/m=%.5f  logq/m=%.4f  delta=%.4f   gain over best pair = %.5f nats/m,  delta gain = %+.4f\n",
     m, btj[1],btj[2],btj[3], log(bt)/m, btq/m, -log(bt)/btq, (log(bp)-log(bt))/m, -log(bt)/btq + log(bp)/bpq),
   printf("m=%d  no admissible triple found\n", m));
}
