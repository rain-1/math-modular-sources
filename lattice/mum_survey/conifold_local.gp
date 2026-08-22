read("cyops.gp");
/* normalised quintic: L = th^4 - z*prod(th+j/5),  MUM z=0, conifold z=1 */
polh = prod(j=1,4, 'th+j/5);
C = vector(5,i,if(i==5,1,0));
for(k=0,4, C[k+1] -= 't*polcoef(polh,k,'th));
P = thetaToD(C);
Ps = shiftD(P,1); RD=recData(Ps);
print("leadPoly = ",factor(leadPoly(RD)),"   d=",RD[1]);
N=24;
mk(u2)={my(imp=vector(N+1,i,0),ini=vector(N+1,i,0));
  imp[1]=1;imp[2]=1;ini[2]=1;imp[3]=1;ini[3]=u2;
  solveSeries(RD,N,ini,imp,[0],0)[1];}
{ my(r0=mk(0),r1=mk(1),G0,G1,K=8,res0,res1,u2);
  G0=logRHS(Ps,r0,N,K); G1=logRHS(Ps,r1,N,K);
  my(imp=vector(N+1,i,0),ini=vector(N+1,i,0)); imp[1]=1;imp[2]=1;imp[3]=1;
  res0=solveSeries(RD,N,ini,imp,vector(#G0[1],i,-G0[1][i]),K)[2];
  res1=solveSeries(RD,N,ini,imp,vector(#G1[1],i,-G1[1][i]),K)[2];
  print("res(u2=0)=",res0,"  res(u2=1)=",res1);
  my(a=res0[1][2], b=res1[1][2]-res0[1][2]); u2=-a/b;
  print("=> canonical u2 = ",u2);
  my(V=mk(u2), GG=logRHS(Ps,V,N,K));
  my(rr=solveSeries(RD,N,ini,imp,vector(#GG[1],i,-GG[1][i]),K));
  print("OBS with canonical u2: ",rr[2]);
  print("V = ",vector(9,i,V[i]));
  print("W = ",vector(9,i,rr[1][i]));
}
quit
