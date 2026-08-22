read("cyops.gp");
/* quintic: L = theta^4 - 5t(5th+1)(5th+2)(5th+3)(5th+4) */
pol = prod(j=1,4, 5*'th+j);
C = vector(5, i, if(i==5, 1, 0));
for(k=0,4, C[k+1] -= 5*'t*polcoef(pol,k,'th));
print("C = ",C);
P = thetaToD(C);
print("D-form degrees: ", vector(#P,i,poldegree(P[i],'t)));
Ps = shiftD(P,0);
RD = recData(Ps);
print("d = ",RD[1], "  leadPoly = ", leadPoly(RD));
N=12;
imp = vector(N+1,i,0); imp[1]=1; ini=vector(N+1,i,0); ini[1]=1;
r = solveSeries(RD,N,ini,imp,[0]);
print("A_n = ", r[1]);
print("check (5n)!/n!^5: ", vector(N+1,i, (5*(i-1))!/(i-1)!^5));
print("OBS=",r[2]," FREE=",r[3]);
quit
