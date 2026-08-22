default(parisizemax,6000000000);
M=400;
A=vector(M+2); A[1]=1; A[2]=5;
{for(n=1,M, A[n+2] = ((2*n+1)*(17*n^2+17*n+5)*A[n+1] - n^3*A[n])/(n+1)^3);}
B=vector(M+2); B[1]=0;
{for(n=1,M+1, my(m=n-1,bm1,bm2,am,am1,am2); bm1=B[n]; bm2=if(n>=2,B[n-1],0); am=A[n+1]; am1=A[n]; am2=if(n>=2,A[n-1],0); B[n+1] = ( (2*n-1)*(17*m^2+17*m+5)*bm1 - m^3*bm2 - 3*n^2*am + (102*m^2+102*m+27)*am1 - 3*m^2*am2 )/n^3 );}
Fs = sum(n=0,M, A[n+1]*t^n) + O(t^(M+1));
Gs = sum(n=0,M, B[n+1]*t^n) + O(t^(M+1));
uu = Gs/Fs;
q  = t*exp(uu);
g  = q/t;
vp(x,p) = if(x==0, 10^6, valuation(x,p));
{for(pi=1,4, my(p=[2,3,5,7][pi], Ftp, Fr, mn, K); K=M-20;
  Ftp = sum(n=0,M\p, A[n+1]*t^(n*p)) + O(t^(M+1));
  Fr = Fs/Ftp; mn=10^6;
  for(i=0,K, mn=min(mn,vp(polcoeff(Fr,i),p)));
  print("D1  p=",p,": min v_p of [t^i] F(t)/F(t^p), i<=",K," : ",mn));}
{for(pi=1,4, my(p=[2,3,5,7][pi], utp, Bf, mn, wi, K); K=M-20;
  utp = sum(n=0,M\p, polcoeff(uu,n)*t^(n*p)) + O(t^(M+1));
  Bf = uu - utp/p; mn=10^6; wi=-1;
  for(i=1,K, my(v=vp(polcoeff(Bf,i),p)); if(v<mn, mn=v; wi=i));
  print("D2  p=",p,": min v_p of [t^i] (u(t)-u(t^p)/p), 1<=i<=",K," : ",mn," at i=",wi));}
{for(pi=1,4, my(p=[2,3,5,7][pi], gtp, h, mn, wi, K); K=M-30;
  gtp = sum(n=0,M\p, polcoeff(g,n)*t^(n*p)) + O(t^(M+1));
  h = gtp/(g^p); mn=10^6; wi=-1;
  for(i=1,K, my(v=vp(polcoeff(h,i),p)); if(v<mn, mn=v; wi=i));
  print("Dwork-h  p=",p,": min v_p of [t^i](g(t^p)/g(t)^p), 1<=i<=",K," : ",mn," at i=",wi));}
quit;
