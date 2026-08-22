read("cyops.gp");
pol = prod(j=1,4, 5*'th+j);
C = vector(5, i, if(i==5, 1, 0));
for(k=0,4, C[k+1] -= 5*'t*polcoef(pol,k,'th));
P = thetaToD(C); tc=1/3125; Ps=shiftD(P,tc); RD=recData(Ps);
N=16;
mk(u2) = { my(imp=vector(N+1,i,0), ini=vector(N+1,i,0));
  imp[1]=1;ini[1]=0; imp[2]=1;ini[2]=1; imp[3]=1;ini[3]=u2;
  solveSeries(RD,N,ini,imp,[0],0)[1]; }
{
for(k=0,3,
  my(u2 = k*1/6, V=mk(u2), GG=logRHS(Ps,V,N,6), G=GG[1], K=GG[2]);
  /* solve L(W) = -G, W free at 0,1,2 -> set to 0 */
  my(imp=vector(N+1,i,0), ini=vector(N+1,i,0));
  imp[1]=1;imp[2]=1;imp[3]=1;
  my(r=solveSeries(RD,N,ini,imp,vector(#G,i,-G[i]),K));
  print("u2=",u2,"  OBS=",r[2]);
);
}
quit
