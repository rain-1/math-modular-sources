read("cyops.gp");
pol = prod(j=1,4, 5*'th+j);
C = vector(5, i, if(i==5, 1, 0));
for(k=0,4, C[k+1] -= 5*'t*polcoef(pol,k,'th));
P = thetaToD(C); tc=1/3125; Ps=shiftD(P,tc); RD=recData(Ps);
N=10;
imp=vector(N+1,i,0); ini=vector(N+1,i,0);
imp[1]=1; ini[1]=0; imp[2]=1; ini[2]=1;
r=solveSeries(RD,N,ini,imp,[0],0);
print("V(u2 free->0) = ",r[1]);
print("OBS=",r[2]," FREE=",r[3]);
/* exponent-0 solution */
imp2=vector(N+1,i,0); ini2=vector(N+1,i,0); imp2[1]=1; ini2[1]=1;
r0=solveSeries(RD,N,ini2,imp2,[0],0);
print("g0 = ",r0[1]); print("OBS=",r0[2]," FREE=",r0[3]);
/* exponent-2 solution */
imp3=vector(N+1,i,0); ini3=vector(N+1,i,0); imp3[1]=1;imp3[2]=1;imp3[3]=1; ini3[3]=1;
r2=solveSeries(RD,N,ini3,imp3,[0],0);
print("g2 = ",r2[1]); print("OBS=",r2[2]," FREE=",r2[3]);
quit
