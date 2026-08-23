/* run_h4.gp -- Task 4: the adelic picture.
   (a) xi2 = lim P_m/Q_m in Q_2 and the exact tail valuation;
   (b) the Nesterenko row converges 2-adically to the SAME xi2;
   (c) the 2-adic linear forms of both rows, and of the output pair (q,p);
   (d) archimedean and 2-adic qualities of the same (q,p).                  */
\p 3000
default(parisize, 6000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/";
RW  = rdrows(concat(DIR,"rows_n200.txt"));
MTOP = 620;
XD = xi2at(MTOP); XI = XD[1]; XPR = XD[2]; XR = XD[3]; ZZ = XD[4];
print("xi2: 2-adic precision (bits) = ", XPR);
print("xi2 valuation = ", valuation(XI));
print("xi2 mod 2^40  = ", lift(XR + O(2^40)));
print("xi2 mod 2^80  = ", lift(XR + O(2^80)));
print("Zudilin tail v_2(xi2 - P_m/Q_m) vs 8m-1-4s_2(m):");
{
for(j=1,6, my(m=40*j, v=valuation(XR - ZZ[2][m+1]/ZZ[1][m+1], 2));
   print("  m=",m,"  v_2=",v,"  8m-1-4s_2(m)=", 8*m-1-4*hammingweight(m),
         "  match=", v==8*m-1-4*hammingweight(m)));
}
print("");
print("k,n,v2ZudForm,v2NestForm,v2Nesttail,v2mixed,v2MOD,logq,archrate,v2qxip,d2,dinf,dad,incone");
KLIST = [22.4, 23.0, 23.9];
{
for(nn=4,200,
 my(rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*nn,i,i)), SS=DD^2,
    vZ = valuation(XX*XR-YY, 2), vN = valuation(VV*XR-UU, 2),
    BN = VV/(4^(7*nn+1)*SS), CN = UU/(4^(7*nn)*SS),
    vNt = valuation(CN - 4*BN*XR, 2),
    vmix = valuation(XX*UU-VV*YY, 2));
 for(ii=1,#KLIST,
  my(kk=KLIST[ii], TT=2^floor(kk*nn), MOD=SS*TT,
     B0=kfull(XX,YY,VV,UU,MOD),
     LZ=(XX*GG-YY)/MOD, LN=(VV*GG-UU)/MOD,
     sZ=sign(LZ), sN=sign(LN), lz=abs(LZ), ln=abs(LN),
     BB=[sZ*B0[1,1], sZ*B0[1,2]; sN*B0[2,1], sN*B0[2,2]],
     gr=gred(BB,lz,ln), Br=gr[1], l1=gr[2],
     cov=abs(matdet(B0))*lz*ln, L0=eight(Br,lz,ln),
     jhi=if(L0==0,8,floor(L0*l1/cov)+2),
     sc=conescan(Br,lz,ln,jhi), cv=sc[2],
     cz=sZ*cv[1], cn=sN*cv[2],
     q=(cz*XX+cn*VV)/MOD, p=(cz*YY+cn*UU)/MOD,
     v1=[Br[1,1],Br[2,1]],
     inc = if((v1[1]>=0&&v1[2]>=0)||(v1[1]<=0&&v1[2]<=0),1,0));
  my(vqp = valuation(q*XR-p, 2), lq=log(abs(q)), la=log(abs(q*GG-p)));
  printf("%.1f,%d,%d,%d,%d,%d,%d,%.5f,%.5f,%d,%.5f,%.5f,%.5f,%d\n",
    kk,nn, vZ, vN, vNt, vmix, valuation(MOD,2), lq/nn, la/nn, vqp,
    vqp*log(2)/lq, -la/lq, (vqp*log(2)-la)/lq, inc)));
}
\q
