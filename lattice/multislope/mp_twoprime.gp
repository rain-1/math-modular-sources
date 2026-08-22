\\ mp_twoprime.gp -- which census rows have TWO slope primes, and what a
\\ two-sampling lattice on such a row would give.  (Fable, 2026-08-22)
default(parisizemax,2000000000);
default(realprecision,60);
NN=600;
zrow(a,b,c,N)={
  my(av=vector(N+2),bv=vector(N+2));
  av[1]=1; av[2]=b; bv[1]=0; bv[2]=1;
  for(n=1,N,
    av[n+2]=((a*n^2+a*n+b)*av[n+1]-c*n^2*av[n])/(n+1)^2;
    bv[n+2]=((a*n^2+a*n+b)*bv[n+1]-c*n^2*bv[n])/(n+1)^2);
  [av,bv]};
slp(av,bv,p,N)={
  my(t=bv[N+1]/av[N+1]-bv[N]/av[N]);
  if(t==0, 0, valuation(t,p)*1./(N-1))};
aa=[7,9,10,11,12,17];
bb=[2,3,3,3,4,6];
cc=[-8,27,9,-1,32,72];
nn=["A","B","C","D","E","F"];
print("row  c    v_p(c) p=2,3,5,7        measured sigma_p at N=",NN);
{
for(i=1,6,
  my(rr=zrow(aa[i],bb[i],cc[i],NN));
  my(vp=vector(4,j,valuation(cc[i],[2,3,5,7][j])));
  my(sp=vector(4,j,slp(rr[1],rr[2],[2,3,5,7][j],NN)));
  print(nn[i],"  ",cc[i],"  ",vp,"   ",sp));
}
print();
print("=== Zagier F: design numbers for a two-sampling lattice on ONE row ===");
{
my(l1=9.0, l2=8.0, k=2, res=3*log(2)+2*log(3));
print("  lambda_1=",l1,"  lambda_2=",l2,"  c=72,  Sum_p sigma_p log p = ",res);
print("  score = log(1/lambda_2)-k = ",log(1/l2)-k);
print("  budget = score + Sum sigma_p log p = ",log(1/l2)-k+res);
for(j=1,5,
  my(al=[0.5,0.75,1.0,1.5,2.0][j], ga=1.0);
  my(G=min(al,ga)*res);
  my(Fv=(k*max(al,ga)+al*log(l2)-ga*log(1/l2)-G)/2);
  my(Hv=Fv+ga*log(l1/l2));
  print("  alpha=",al,"  G=",G,"  F=",Fv,"  H=",Hv,"  delta=",1-Fv/Hv));
}
\q
