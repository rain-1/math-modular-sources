read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/ode_fit.gp");
NT=80; r=zudrow(NT); u=r[1]; v=r[2];
U = vector(NT+1, i, 16^(i-1)*u[i]);
y = sum(n=0,NT-2, U[n+1]*x^n);
{ for(rr=2,3, for(dm=2,10,
   my(d=matsize(odefit(y,rr,dm)[1])[2]);
   print("order ",rr," deg ",dm,": ker dim ",d))); }
PR=400; N=60;
xi = v[N+1]/u[N+1]; xim = v[N]/u[N];
print("N=",N-1,"  v_2(xi_N - xi_{N-1}) = ", valuation(xi-xim,2));
z2 = Lp(2,triv,2,PR);
print("v_2(xi - zeta_2(2))   = ", valuation(xi - z2,2));
print("v_2(xi - zeta_2(2)/2) = ", valuation(xi - z2/2,2));
print("v_2(zeta_2(2)) = ", valuation(z2,2));
quit;
