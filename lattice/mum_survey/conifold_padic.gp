read("cyops.gp");
polh = prod(j=1,4,'th+j/5);
CC = vector(5,i,if(i==5,1,0));
for(k=0,4, CC[k+1] -= 't*polcoef(polh,k,'th));
P = thetaToD(CC); Ps = shiftD(P,1); RD = recData(Ps);
NC = 260;
imp=vector(NC+1,i,0); ini=vector(NC+1,i,0);
imp[1]=1;imp[2]=1;ini[2]=1;imp[3]=1;ini[3]=-7/10;
V = solveSeries(RD,NC,ini,imp,[0],0)[1];
GG = logRHS(Ps,V,NC,8);
imp2=vector(NC+1,i,0); ini2=vector(NC+1,i,0); imp2[1]=1;imp2[2]=1;imp2[3]=1;
rr = solveSeries(RD,NC,ini2,imp2,vector(#GG[1],i,-GG[1][i]),GG[2]);
W = rr[1];
print("OBS = ",rr[2]);
print("v_p(v_n) rates:");
{ for(pi=1,4, my(p=[2,3,5,7][pi]);
   print("  p=",p,": ",vector(6,i,my(n=40*i); valuation(V[n+1],p)));
  );
}
print("v_p(w_n/v_n - w_{n-1}/v_{n-1}):");
{ for(pi=1,4, my(p=[2,3,5,7][pi]);
   print("  p=",p,": ",vector(8,i,my(n=30*i); valuation(W[n+1]/V[n+1]-W[n]/V[n],p)));
  );
}
default(realprecision,40);
print("w_n/v_n (arch): ", vector(6,i,my(n=40*i); 1.0*W[n+1]/V[n+1]));
quit
