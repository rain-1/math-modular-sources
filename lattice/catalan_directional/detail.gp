\p 60
{
detail(n) = my(D=lcm(vector(6*n,i,i)), zr=zudrow(n), nr=nestrow(n),
   X=zr[1],Y=zr[2],V=nr[1],U=nr[2], Sw=1);
 forprime(l=2,6*n, if((l>4*n&&l<6*n)||(l>2*n&&l<3*n), Sw*=l^2));
 printf("--- n=%d  (D^2 mass=%.4f/n, window Sw mass=%.4f/n)\n", n, 2*log(D)/n, log(Sw)/n);
 printf("  l  vX vY vV vU | s(D^2) eOpt fZ(D^2) net | s(win) fZ(win) net_win\n");
 forprime(l=3,6*n,
  my(vX=valuation(X,l),vY=valuation(Y,l),vV=valuation(V,l),vU=valuation(U,l),
     s=2*valuation(D,l), sw=valuation(Sw,l));
  my(eO=max(0,min(vV,vU)-min(vX,vY)));
  if(eO>0 || vX!=vY || vV!=vU,
   printf("%4d  %2d %2d %2d %2d |  %2d    %2d    %2d    %2d |  %2d    %2d    %2d\n",
     l,vX,vY,vV,vU,s,eO,max(0,min(vU,s)-vY),max(0,eO-max(0,min(vU,s)-vY)),
     sw,max(0,min(vU,sw)-vY),max(0,eO-max(0,min(vU,sw)-vY)))));
}
detail(20);
