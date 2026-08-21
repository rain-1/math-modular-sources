read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/ode_fit.gp");
NT=70; r=zudrow(NT); u=r[1];
U = vector(NT+1, i, 16^(i-1)*u[i]);
print("U_n = 16^n u_n : ", vector(10,i,U[i]));
print("integral? ", vector(NT+1,i,denominator(U[i]))==vector(NT+1,i,1));
y = sum(n=0,NT-2, U[n+1]*x^n);
{ for(rr=2,4, for(dm=2,6,
   my(d=matsize(odefit(y,rr,dm)[1])[2]);
   if(d>0, print("order ",rr," deg ",dm,": ker dim ",d)))); }
{ my(z=odefit(y,2,4)); if(matsize(z[1])[2]>0, print("order2 deg4 relation:"); showrel(z[1],z[2])); }
{ my(z=odefit(y,2,5)); if(matsize(z[1])[2]>0, print("order2 deg5 relation:"); showrel(z[1],z[2])); }
{ my(z=odefit(y,3,4)); if(matsize(z[1])[2]>0, print("order3 deg4 relation:"); showrel(z[1],z[2])); }
quit;
