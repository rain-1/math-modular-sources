/* run_h3.gp: where G enters.  For truncations G_d of G at d decimal digits,
   recompute the ORIENTATION, the Hermite data, the whole continued fraction
   of h12/h11, and the balance index.  Reports for each (n,d): whether the
   sign sZ, the partial-quotient string, and the index agree with the truth. */
\p 3000
default(parisize, 4000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW  = rdrows(concat(DIR,"rows_all.txt"));
NLIST = [8, 16, 24, 40, 60, 80, 100, 120];
KK = 23.9;
print("n,d,dstar,signok,cfok,idx,idxtrue,dcf");
{
for(jj=1,#NLIST,
 my(nn=NLIST[jj], rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*nn,i,i)), SS=DD^2,
    hd=hdat(XX,YY,VV,UU,SS,KK,nn,0), h11=hd[1], h12=hd[2], h22=hd[3], MOD=hd[4],
    ld=ladder(h11,h12,h22),
    t0=((VV*GG-UU)/abs(XX*GG-YY))^2, b0=balance(ld[1],ld[2],t0),
    cf0=ld[3], ds=6.2675*nn);
 my(dl=List());
 forstep(d=floor(ds)-30, floor(ds)+30, 6, listput(dl,d));
 listput(dl, floor(ds)+120); listput(dl, floor(ds)+400);
 for(m=1,#dl,
  my(d=dl[m], Gd=floor(GG*10^d)/10^d*1.0,
     lz0=(XX*Gd-YY)/MOD, ln0=(VV*Gd-UU)/MOD,
     sZd=sign(lz0), sNd=sign(ln0));
  my(B0=kfull(XX,YY,VV,UU,MOD),
     BBd=[sZd*B0[1,1], sZd*B0[1,2]; sNd*B0[2,1], sNd*B0[2,2]],
     HHd=mathnf(BBd), ldd=ladder(HHd[1,1],HHd[1,2],HHd[2,2]),
     td=(abs(ln0)/abs(lz0))^2, bd=balance(ldd[1],ldd[2],td),
     cfd=ldd[3], same=(cfd==cf0), dc=0);
  if(!same, for(i=1,min(#cfd,#cf0), if(cfd[i]!=cf0[i], dc=i; break)));
  printf("%d,%d,%.1f,%d,%d,%d,%d,%d\n", nn, d, ds,
     sZd==hd[6], if(same,1,0), bd[1], b0[1], dc)));
}
\q
