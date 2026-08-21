read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/ode_fit.gp");
NT=50;
A=vector(NT+2); A[1]=1; A[2]=6;
for(n=1,NT, A[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*A[n+1] - 12*n*(16*n^2-1)*A[n])/(n+1)^3);
y = sum(n=0,NT, A[n+1]*x^n);
{
for(rr=2,3, for(dm=2,4,
  my(z=odefit(y,rr,dm));
  print("s18: order ",rr," deg ",dm," kernel dim ", matsize(z[1])[2]);
  if(matsize(z[1])[2]>0 && rr==2, showrel(z[1],z[2]))));
}
quit;
