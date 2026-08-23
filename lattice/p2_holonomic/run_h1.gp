\p 200
default(parisize, 2000000000);
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW  = rdrows(concat(DIR,"rows_all.txt"));
KLIST = [22.4, 23.0, 23.9];
print("k,n,h11ok,h22ok,okY,okX,lcmok,logh11,logh22,logmY,logmX,v2h11,logoddh11,logmZ,v2mZ,logg0,v2gMODY,logh12");
{
for(nn=4,120,
 my(rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*nn,i,i)), SS=DD^2);
 for(ii=1,#KLIST,
  my(kk=KLIST[ii], hd=hdat(XX,YY,VV,UU,SS,kk,nn,0),
     h11=hd[1], h12=hd[2], h22=hd[3], MOD=hd[4], idx=hd[5], sZ=hd[6],
     pr=hpred(XX,YY,VV,UU,MOD),
     ck=hchk12(XX,YY,VV,UU,MOD,h11,h12,h22,sZ),
     mZ=gcd(h11,h12));
  printf("%.1f,%d,%d,%d,%d,%d,%d,%.5f,%.5f,%.5f,%.5f,%d,%.5f,%.5f,%d,%.5f,%d,%.5f\n",
    kk,nn, h11==pr[1], h22==pr[2], ck[1],ck[2],ck[5],
    log(h11)/nn, log(h22)/nn, log(ck[3])/nn, log(ck[4])/nn,
    valuation(h11,2), log(h11/2^valuation(h11,2))/nn,
    log(mZ)/nn, valuation(mZ,2), log(pr[3])/nn,
    valuation(gcd(MOD,YY),2), log(h12)/nn)));
}
\q
