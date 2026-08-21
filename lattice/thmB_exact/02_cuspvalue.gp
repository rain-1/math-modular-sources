/* 02_cuspvalue.gp -- the archimedean endpoint criterion and the cusp value of
   the Eichler integral.

   Theta(iy) = (polynomial of degree w in y) + O(exp(-c/y)) whose constant term
   is L(Phi,w+1), PROVIDED Phi vanishes at the cusp 0, i.e. delta(phi)P(w+2)=0.
   Richardson in y removes the higher terms of the polynomial.              */
read("common.gp");
default(realprecision, 70);
Y0 = 0.004;
MM = 26000;
Th(cv,y) = my(u=exp(-2*Pi*y), M=min(MM,ceil(70*log(10)/(2*Pi*y))), s=0.); forstep(m=M,1,-1, s=(s+cv[m])*u); s;
print("row      d(phi)P(w+2)   Theta(0) - L(Phi,w+1)");
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], w=r-1);
  my(cv=vector(MM,m,cm(nm,m)*1.0/m^(w+1)));
  my(f0=Th(cv,Y0), f1=Th(cv,Y0/2), f2=Th(cv,Y0/4), L=Ltarget(nm), ex);
  ex = if(w==1, 2*f1-f0, (8/3)*f2-2*f1+(1/3)*f0);
  print(nm,"\t",endpoint(nm,r),"\t",ex-L, if(endpoint(nm,r)==0,"","   <- diverges: no finite cusp value"));
);
}
quit;
