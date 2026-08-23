/* run_h7.gp -- does the 2-adic modulus control the continued fraction?
   Fix n and slide the 2-adic exponent E = v_2(T_n) by one at a time.  Each
   step multiplies the modulus by 2 and h11 by (at most) 2.  If there were a
   2-adic law behind the CF of h12/h11, one step should act predictably.     */
\p 3000
default(parisize, 6000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/";
RW  = rdrows(concat(DIR,"rows_n200.txt"));
NL = [20, 40, 60, 80, 100, 120];
print("n,E,logh11,len,idx,parity,a1,a2,a3,cfagree,pqnext");
{
for(jj=1,#NL,
 my(nn=NL[jj], rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*nn,i,i)), SS=DD^2, E0=floor(22.4*nn), prev=0);
 for(E=E0, E0+60,
  my(MOD=SS*2^E, B0=kfull(XX,YY,VV,UU,MOD), sZ=(-1)^nn,
     BB=[sZ*B0[1,1], sZ*B0[1,2]; B0[2,1], B0[2,2]], HH=mathnf(BB),
     h11=HH[1,1], h12=HH[1,2], h22=HH[2,2],
     LZ=abs(XX*GG-YY)/MOD, LN=(VV*GG-UU)/MOD, t=(LN/LZ)^2,
     ld=ladder(h11,h12,h22), bal=balance(ld[1],ld[2],t), cf=ld[3],
     ag=0);
  if(prev!=0, for(i=1,min(#cf,#prev), if(cf[i]!=prev[i], ag=i-1; break));
              if(ag==0, ag=min(#cf,#prev)));
  printf("%d,%d,%.5f,%d,%d,%d,%d,%d,%d,%d,%d\n", nn, E, log(h11)/nn, #cf,
     bal[1], bal[1]%2, if(#cf>1,cf[2],0), if(#cf>2,cf[3],0), if(#cf>3,cf[4],0),
     ag, if(bal[1]+2<=#cf, cf[bal[1]+2], 0));
  prev = cf));
}
\q
