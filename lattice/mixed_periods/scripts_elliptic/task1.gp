default(parisize,"2G");
N=300;
lc=vector(N+3); lc[1]=1; for(i=2,N+3, lc[i]=lcm(lc[i-1],i));
LC(n)=if(n<1,1,lc[n]);

excess(a,D)=my(d=denominator(a)); d/gcd(d,D);

analyse(name,av)=
{
  my(PR=N+3, H, KB, KD, KL, HB, HD, HL, T, mx, mxn, nonint, c);
  H = 1+O(x^PR);
  for(j=1,#av, H = H*(1-av[j]*x+O(x^PR))^(-1/2));
  \\ integrality of H
  nonint=0; for(n=0,N, if(denominator(polcoeff(H,n))!=1, nonint++));
  print(name, "  H integral for n<=",N,": ", if(nonint==0,"YES", Str("NO, ",nonint," bad")));
  KB = 1/(1-x+O(x^PR));
  KL = log(1-x+O(x^PR));
  KD = KL*KB;
  HB = H*intformal(H*KB);
  HD = H*intformal(H*KD);
  HL = H*intformal(H*KL);
  \\ H_B vs [1..n]
  mx=1; mxn=-1; for(n=1,N, c=excess(polcoeff(HB,n), LC(n)); if(c>mx, mx=c; mxn=n));
  print("   H_B  max excess over lcm(1..n)            = ",mx, if(mxn>=0, Str("  at n=",mxn), ""));
  mx=1; mxn=-1; for(n=1,N, c=excess(polcoeff(HD,n), LC(n)*LC(n\2)); if(c>mx, mx=c; mxn=n));
  print("   H_D  max excess over lcm(1..n)lcm(1..n/2) = ",mx, if(mxn>=0, Str("  at n=",mxn), ""));
  mx=1; mxn=-1; for(n=1,N, c=excess(polcoeff(HD,n), LC(n)); if(c>mx, mx=c; mxn=n));
  print("        (max excess over lcm(1..n) alone     = ",mx, if(mxn>=0, Str("  at n=",mxn), ""),")");
  mx=1; mxn=-1; for(n=1,N, c=excess(polcoeff(HL,n), LC(n)*LC(n\2)); if(c>mx, mx=c; mxn=n));
  print("   H_L  max excess over lcm(1..n)lcm(1..n/2) = ",mx, if(mxn>=0, Str("  at n=",mxn), ""));
  mx=1; mxn=-1; for(n=1,N, c=excess(polcoeff(HL,n), LC(n)); if(c>mx, mx=c; mxn=n));
  print("        (max excess over lcm(1..n) alone     = ",mx, if(mxn>=0, Str("  at n=",mxn), ""),")");
  print("   first coeffs H: ", vector(8,i,polcoeff(H,i-1)));
  print("");
}

for(k=1,7, my(mv=[1,2,3,4,5,6,12], m=mv[k]); analyse(Str("E_",m), [1,9,4*m]));
for(m=1,3, analyse(Str("F_",m), [1,25,4*m]));
analyse("G", [4,8,12]);
quit;
