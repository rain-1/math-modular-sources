/* run_h8.gp -- is the E-persistence of the parity a property of the Catalan
   rows, or of index-2 sublattice geometry?  Same slide, with a control in
   which (X,Y,V,U) are replaced by random integers of the same 2-adic and
   archimedean sizes (RND=1), and with the lag structure recorded.          */
\p 2000
default(parisize, 6000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/";
RW  = rdrows(concat(DIR,"rows_n200.txt"));
NL = [20,30,40,50,60,70,80,90,100,110,120];
print("tag,n,E,idx,parity,survive,logl1,logl2,skew");
{
slide(nn, XX,YY,VV,UU, SS, tag, ESPAN) =
 my(E0=floor(22.4*nn), prevv=0);
 for(E=E0, E0+ESPAN,
  my(MOD=SS*2^E, B0=kfull(XX,YY,VV,UU,MOD),
     LZ=(XX*GG-YY)/MOD, LN=(VV*GG-UU)/MOD, sZ=sign(LZ), sN=sign(LN),
     lz=abs(LZ), ln=abs(LN),
     BB=[sZ*B0[1,1], sZ*B0[1,2]; sN*B0[2,1], sN*B0[2,2]],
     HH=mathnf(BB), h11=HH[1,1], h12=HH[1,2], h22=HH[2,2],
     gr=gred(BB,lz,ln), l1=gr[2], l2=gr[3],
     ld=ladder(h11,h12,h22), bal=balance(ld[1],ld[2],(ln/lz)^2),
     vv=[ld[1][bal[1]+1], ld[2][bal[1]+1]],
     sv=if(prevv==0, -1, if(vv==prevv, 1, 0)));
  printf("%s,%d,%d,%d,%d,%d,%.5f,%.5f,%.4f\n", tag, nn, E, bal[1], bal[1]%2, sv,
     log(l1)/nn, log(l2)/nn, l2/l1);
  prevv = vv);
}
{
for(jj=1,#NL,
 my(nn=NL[jj], rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*nn,i,i)), SS=DD^2);
 slide(nn, XX,YY,VV,UU, SS, "cat", 120));
}
/* control: random rows with the same sizes and the same 2-adic cross
   divisibility switched OFF (so only sizes match)                         */
setrand(20260823);
{
for(jj=1,#NL,
 my(nn=NL[jj], rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*nn,i,i)), SS=DD^2,
    rX=2*random(XX)+1, rV=2*random(VV)+1,
    rY=round(rX*GG) - random(floor(abs(XX*GG-YY))),
    rU=round(rV*GG) - random(floor(abs(VV*GG-UU))));
 slide(nn, rX,rY,rV,rU, SS, "rnd", 120));
}
\q
