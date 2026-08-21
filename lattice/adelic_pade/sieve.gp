/* Adelic sieve on the two-row lattice: 2-adic margin of selected vectors.
   margin(c) := v2(G2 - p/q).  Kernel vectors (rational G) have margin O(1);
   the nonvanishing criterion is margin -> +oo. */
s2(x)=hammingweight(x);
{
klat2(X,Y,V,U,M) = my(K=matkerint([X,V,M,0,0,0; Y,U,0,M,0,0]));
  mathnf(matrix(2,#K[1,],i,j,K[i,j]));
}
{
redu(B,lz,ln) = my(sc=10^(400+ceil(log(1/min(lz,ln))/log(10))));
  B*qflll([round(B[1,1]*lz*sc),round(B[1,2]*lz*sc); round(B[2,1]*ln*sc),round(B[2,2]*ln*sc)]);
}
{
run(nlist,k,R,GG,G2) = my();
 printf("  n |  logL1/n | margin(l1) | margin(l2) | 24n-1-4s2(3n) | 28n-.. | marg(cone) | minmarg over ball | #short with marg<24n-1\n");
 for(t=1,#nlist,
  my(n=nlist[t], D=lcm(vector(6*n,i,i)), S=D^2, T=2^floor(k*n), M=S*T,
     zr=zudrow(n), nr=nestrow(n), X=zr[1],Y=zr[2],V=nr[1],U=nr[2]);
  my(B0=klat2(X,Y,V,U,M));
  my(LZ=(X*GG-Y)/M, LN=(V*GG-U)/M, sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN));
  my(B=[sZ*B0[1,1],sZ*B0[1,2]; sN*B0[2,1],sN*B0[2,2]], Br=redu(B,lz,ln));
  my(mrg(cz,cn)= my(q=(cz*X+cn*V)/M,p=(cz*Y+cn*U)/M); if(q==0,-1,v2(G2-p/q)));
  my(v1=[sZ*Br[1,1],sN*Br[2,1]], v2v=[sZ*Br[1,2],sN*Br[2,2]]);
  my(m1=mrg(v1[1],v1[2]), m2=mrg(v2v[1],v2v[2]));
  my(qZ=24*n-1-4*s2(3*n), qN=28*n-1-2*s2(n)-2*s2(3*n));
  /* scan a ball of short vectors, record margins */
  my(mm=10^9, cnt=0, tot=0, bestcone=0, mcone=-1);
  for(i=-R,R, for(jj=-R,R, my(cz=i*Br[1,1]+jj*Br[1,2], cn=i*Br[2,1]+jj*Br[2,2]);
    if(cz||cn, my(q=(cz*sZ*X+cn*sN*V)/M);
      if(q!=0, my(mg=mrg(sZ*cz,sN*cn)); tot++; if(mg<mm,mm=mg); if(mg<qZ,cnt++);
        if(cz>=0&&cn>=0, my(f=cz*lz+cn*ln); if(bestcone==0||f<bestcone, bestcone=f; mcone=mg))))));
  printf("%4d | %8.4f | %10d | %10d | %13d | %6d | %10d | %8d | %d/%d\n",
    n, log(sqrt((v1[1]*lz)^2+(v1[2]*ln)^2))/n, m1, m2, qZ, qN, mcone, mm, cnt, tot);
 );
}
G2=G2rat(140)[1];
print("G2 precision 2^",G2rat(140)[2]);
run([4,6,8,10,12,14,16],22.4,25,Catalan,G2);
\q
