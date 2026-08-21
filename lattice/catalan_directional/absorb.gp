\p 60
/* Checks: (i) v_l(Y_n)=0 for all odd l<=6n; (ii) v_l(U_n)<=v_l(S_n)=2v_l(D_{6n});
   (iii) hence eOpt = fZ (total absorption) under S=D^2.
   Also measures net directional mass under the literal-window modulus S_w (3.6). */
{
for(idx=1,14, my(n=[5,6,8,10,12,15,18,20,24,28,32,36,44,60][idx]);
 my(D=lcm(vector(6*n,i,i)), zr=zudrow(n), nr=nestrow(n),
    X=zr[1],Y=zr[2],V=nr[1],U=nr[2], Sw=1, badY=0, badU=0, netD=0.0, netW=0.0, opt=0.0);
 forprime(l=2,6*n, if((l>4*n&&l<6*n)||(l>2*n&&l<3*n), Sw*=l^2));
 forprime(l=3,6*n,
  my(vX=valuation(X,l),vY=valuation(Y,l),vV=valuation(V,l),vU=valuation(U,l),
     s=2*valuation(D,l), sw=valuation(Sw,l), L=log(l));
  if(vY!=0, badY++); if(vU>s, badU++);
  my(eO=max(0,min(vV,vU)-min(vX,vY)));
  opt += eO*L;
  netD += max(0,eO-max(0,min(vU,s)-vY))*L;
  netW += max(0,eO-max(0,min(vU,sw)-vY))*L);
 printf("n=%2d  #{l odd: v_l(Y)!=0}=%d  #{v_l(U)>v_l(S)}=%d | optmass/n=%.5f  net(S=D^2)/n=%.5f  net(S=window)/n=%.5f  [logSw/n=%.4f]\n",
   n,badY,badU,opt/n,netD/n,netW/n,log(Sw)/n));
}
