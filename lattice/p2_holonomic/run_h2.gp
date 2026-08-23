/* run_h2.gp: the balance index i(n), its exact stability interval in the
   metric ratio, and the number of digits of G that decide it.            */
\p 3000
default(parisize, 4000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW  = rdrows(concat(DIR,"rows_all.txt"));
KLIST = [22.4, 23.0, 23.9];
print("k,n,idx,len,parity,logt,logtlo,logthi,marglo,marghi,marg,dstar,logh11,logh22,pqnext,pqat");
{
for(nn=4,120,
 my(rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*nn,i,i)), SS=DD^2);
 for(ii=1,#KLIST,
  my(kk=KLIST[ii], hd=hdat(XX,YY,VV,UU,SS,kk,nn,0),
     h11=hd[1], h12=hd[2], h22=hd[3], MOD=hd[4],
     LZ=abs(XX*GG-YY)/MOD, LN=(VV*GG-UU)/MOD, t=(LN/LZ)^2,
     ld=ladder(h11,h12,h22), bal=balance(ld[1],ld[2],t),
     i0=bal[1], tlo=bal[2], thi=bal[3],
     ml=(log(t)-log(tlo))/2, mh=if(thi==oo, 1e9, (log(thi)-log(t))/2),
     mg=min(ml,mh),
     /* delta of G that moves log|XG-Y| by mg :  delta*X/|XG-Y| = mg      */
     dstar = -log(mg*abs(XX*GG-YY)/XX)/log(10));
  printf("%.1f,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.5f,%.5f,%d,%d\n",
    kk,nn,i0,#ld[1],i0%2, log(t), log(tlo), if(thi==oo,0,log(thi)),
    ml, mh, mg, dstar, log(h11)/nn, log(h22)/nn,
    if(i0+2<=#ld[3], ld[3][i0+2], 0), if(i0+1<=#ld[3], ld[3][i0+1], 0))));
}
\q
