default(realprecision,40);
N=2500; a=17;b=6;c=72;
A0=1;A1=b;B0=0;B1=1;
for(n=1,N-1,A2=((a*n^2+a*n+b)*A1-c*n^2*A0)/(n+1)^2;B2=((a*n^2+a*n+b)*B1-c*n^2*B0)/(n+1)^2;A0=A1;A1=A2;B0=B1;B1=B2);
l=B1*1./A1; print("F lim=",l,"  err~",(8/9.)^N);
G=sumalt(k=0,(-1)^k/(2*k+1)^2); L3=lfun(lfuncreate(-3),2);
print(lindep([l,1,zeta(2),G,L3,log(2),log(3),log(2)^2,log(3)^2,log(2)*log(3),Pi*log(2),Pi*log(3)])~);
\q
