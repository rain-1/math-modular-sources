read("cyops.gp");
pol = prod(j=1,4, 5*'th+j);
C = vector(5, i, if(i==5, 1, 0));
for(k=0,4, C[k+1] -= 5*'t*polcoef(pol,k,'th));
P = thetaToD(C);
tc = 1/3125;
Ps = shiftD(P,tc);
RD = recData(Ps);
print("d = ",RD[1]);
print("leadPoly(M) = ", factor(leadPoly(RD)));
N=10;
/* exponent-1 solution: u_0=0, u_1=1 imposed */
imp=vector(N+1,i,0); ini=vector(N+1,i,0);
imp[1]=1; ini[1]=0; imp[2]=1; ini[2]=1;
r=solveSeries(RD,N,ini,imp,[0]);
print("V = ",r[1]);
print("OBS=",r[2]," FREE=",r[3]);
quit
