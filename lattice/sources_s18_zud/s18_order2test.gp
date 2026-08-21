read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/ode_fit.gp");
NT=90;
A=vector(NT+2); A[1]=1; A[2]=6;
for(n=1,NT, A[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*A[n+1] - 12*n*(16*n^2-1)*A[n])/(n+1)^3);
y = sum(n=0,NT, A[n+1]*x^n);
{ for(dm=2,10, print("order 2, deg ",dm,": ker dim ", matsize(odefit(y,2,dm)[1])[2])); }
{ my(z=odefit(y,3,2)); print("order 3 deg 2 relation:"); showrel(z[1],z[2]); }
quit;
