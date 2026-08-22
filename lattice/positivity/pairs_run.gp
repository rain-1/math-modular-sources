/* driver: all pairs (j1<j2) at a fixed m, both with and without the 2-adic bridge */
{
runm(m) =
 my(NN=m+1, AA=vector(NN), BB=vector(NN), MM=vector(NN));
 for(j=0,m, my(r=mom(m,j)); AA[j+1]=r[1]; BB[j+1]=r[2]; MM[j+1]=r[1]*Catalan-r[2]);
 for(j1=0,m-1, for(j2=j1+1,m,
   my(s=pairstat(AA[j1+1],BB[j1+1],AA[j2+1],BB[j2+1],MM[j1+1],MM[j2+1],1,"br"));
   if(s[1]!=0,
    printf("%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f\n",
      m,j1,j2, log(s[1])/m, log(s[2])/m, s[3], s[4]/m, s[5], s[6]*log(2)/m, s[7]/m))));
}
