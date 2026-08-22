default(realprecision,60);
A=readvec("lattice/root_rows/rows_zeta5_L16.txt");
f(t)=sum(n=0,398,A[n+1]*t^n);
{my(a=0.0212,b=0.0213); for(i=1,300,my(c=(a+b)/2); if(f(a)*f(c)<=0,b=c,a=c));
 my(z=(a+b)/2); print("zero t0 = ",z); print("1/t0 = ",1/z);
 print("algdep(1/t0,2..6): ",vector(5,j,algdep(1/z,j+1)));}
\q
