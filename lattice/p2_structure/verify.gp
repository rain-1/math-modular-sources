/* driver: exact identities used in P2_STRUCTURE.md.
   cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
       lattice/p2_structure/verify.gp > run.gp && gp -q run.gp             */
\p 200
default(parisize, 2000000000);
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW = rdrows(concat(DIR,"rows_all.txt"));
{
chk(NHI) = my(b1=0,b2=0,b3=0,b4=0,b5=0,acc=0.,cnt=0);
 for(n=4,NHI,
  my(rw=mapget(RW,n), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
     SS=lcm(vector(6*n,i,i))^2, g0=gcd(gcd(XX,YY),gcd(VV,UU)));
  if(XX%SS || VV%SS, b1++);
  foreach([22.4,23.0,23.9], kk,
   my(TT=2^floor(kk*n), MOD=SS*TT, B0=kfull(XX,YY,VV,UU,MOD), idx=abs(matdet(B0)));
   cnt++;
   if(MOD%idx, b2++);
   if(MOD/idx != g0, b3++);
   /* content of K_n and its Zudilin-directional content m_Z = gcd(h11,h12) */
   my(HH=mathnf(B0), mZ=gcd(HH[1,1],HH[1,2]));
   if(gcd(gcd(HH[1,1],HH[1,2]),HH[2,2]) != 1, b4++);
   if(mZ<=1, b5++);
   acc += log(1.*mZ)/n));
 printf("instances %d: S nondividing rows %d | idx not dividing MOD %d | MOD/idx != gcd(X,Y,V,U) %d | content(K_n) != 1 %d | m_Z = 1 %d | mean log(m_Z)/n %.4f\n",
   cnt, b1, b2, b3, b4, b5, acc/cnt);
}
chk(120);

/* the kernel output identity  (q,p) = -(h_n/T_n)(b,a)  for a surrogate a/b   */
{
qpchk(NHI) = my(bad=0, cnt=0, GG=Catalan);
 for(e=2,5, my(GS=bestappr(GG,10^e), aa=numerator(GS), bb=denominator(GS));
  for(n=6,NHI, my(rw=mapget(RW,n), XX=rw[1],YY=rw[2],VV=rw[3],UU=rw[4],
     SS=lcm(vector(6*n,i,i))^2, hh=(XX*UU-YY*VV)/SS);
   foreach([22.4,23.9], kk,
    my(TT=2^floor(kk*n), MOD=SS*TT, ck=[aa*VV-bb*UU, -(aa*XX-bb*YY)]);
    my(q=(ck[1]*XX+ck[2]*VV)/MOD, p=(ck[1]*YY+ck[2]*UU)/MOD);
    cnt++;
    if(q != -bb*hh/TT || p != -aa*hh/TT || denominator(hh/TT)!=1, bad++))));
 printf("kernel output (q,p) = -(h_n/T_n)(b,a): %d instances, %d failures\n", cnt, bad);
}
qpchk(30);
\q
