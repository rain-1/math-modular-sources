/* lattice/p2_holonomic/verify.gp -- the exact identities of P2_HOLONOMIC.md Prop.1.
   Prepend rows_pos.gp, p2core.gp, hcore.gp.  Output: 0 failures in 591 instances. */
\p 50
default(parisize, 4000000000);
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/";
RW  = rdrows(concat(DIR,"rows_n200.txt"));
KLIST=[22.4,23.0,23.9]; BAD1=0; BAD2=0; BAD3=0; TOT=0;
{
for(nn=4,200,
 my(rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    DD=lcm(vector(6*nn,i,i)), SS=DD^2);
 for(ii=1,#KLIST,
  my(kk=KLIST[ii], hd=hdat(XX,YY,VV,UU,SS,kk,nn,0),
     h11=hd[1],h12=hd[2],h22=hd[3],MOD=hd[4],sZ=hd[6],
     gY=gcd(MOD,YY), mY=MOD/gY, g=gcd(MOD,gcd(XX,YY)));
  TOT++;
  if(mY!=h11, BAD1++);
  if(gY!=g, BAD2++);
  /* explicit modular-inverse form: h12 = -sZ*h22*(UU/gY)*(YY/gY)^-1 mod h11 */
  if((h22*UU)%gY!=0, BAD3++,
    my(pred = lift(Mod(-sZ*(h22*UU/gY)*lift(Mod(YY/gY,h11)^-1), h11)));
    if(pred!=h12, BAD3++))));
}
print("instances ", TOT);
print("M/gcd(M,Y) == h11 failures:            ", BAD1);
print("gcd(M,Y) == gcd(M,X,Y) failures:       ", BAD2);
print("h12 = -sZ h22 U Y^-1 mod h11 failures: ", BAD3);
\q
