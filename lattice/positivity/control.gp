/* lattice/positivity/control.gp --- the rational control of CATALAN_AUDIT sec.4(a),
   applied to every new table in POSITIVITY_PROGRAM.md.  Replace G by the rational
   G* = bestappr(G, 10^320) and rerun.  Any column that is bit-identical is, by
   construction, not evidence about G.  Prepend rows_pos.gp, cone80.gp, pairs.gp. */
{
ctrl_pairs(m, GS) =
 my(NN=m+1, AA=vector(NN), BB=vector(NN), M1=vector(NN), M2=vector(NN));
 for(j=0,m, my(r=mom(m,j)); AA[j+1]=r[1]; BB[j+1]=r[2];
            M1[j+1]=r[1]*Catalan-r[2]; M2[j+1]=r[1]*GS-r[2]);
 my(bp1=0,bj1=0,bp2=0,bj2=0);
 for(j1=0,m-1, for(j2=j1+1,m,
   my(s=pairstat(AA[j1+1],BB[j1+1],AA[j2+1],BB[j2+1],M1[j1+1],M1[j2+1],1,""));
   if(s[5]==1 && (bp1==0 || s[1]<bp1), bp1=s[1]; bj1=[j1,j2]);
   my(t=pairstat(AA[j1+1],BB[j1+1],AA[j2+1],BB[j2+1],M2[j1+1],M2[j2+1],1,""));
   if(t[1]!=0 && (bp2==0 || t[1]<bp2), bp2=t[1]; bj2=[j1,j2])));
 printf("m=%2d  true G: best (%d,%d) rate/m=%.6f    |   G*=bestappr: best (%d,%d) rate/m=%.6f   %s\n",
   m, bj1[1],bj1[2], log(bp1)/m, bj2[1],bj2[2], log(bp2)/m,
   if(abs(bp1/bp2-1) < 1e-100, "IDENTICAL to >100 digits", "differ"));
}
