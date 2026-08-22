/* Does the level-16 host also carry fold-regular level-32 sources (V8 E)? */
read("00_setup.gp");
G = Catalan;
H16 = mkhost(x16,F16); Av = Aof(H16); NMAX = NT-2;
BB = vector(4, j, Bof(H16, mkPhi(Sin,[[2^(j-1),1]])));   /* E, V2E, V4E, V8E */
rate(S,n) = abs(polcoeff(S,n))^(1./n);
P2v(c) = sum(j=1,#c, c[j]*(2^(j-1))^(-2));
{ TS = [[1,0,0,0],[0,1,-8,0],[0,0,1,-8],[0,1,0,0],[0,0,0,1],[0,0,1,0],[1,4,-32,0],[0,1,4,-32]]; }
{for(i=1,#TS, my(c=TS[i], Bv=sum(j=1,4,c[j]*BB[j]), xi=-P2v(c)*G/2, Cv=Bv-xi*Av);
  print("  c=",c,"  b/a(",NMAX,")=",polcoeff(Bv,NMAX)/polcoeff(Av,NMAX)*1.,
        "  pred=",xi,"  |c_n|^1/n=",rate(Cv,NMAX),"  v2(b_n)/n=",valuation(polcoeff(Bv,NMAX),2)/(NMAX*1.)));}
