/* eichler.gp -- Theta(iy) -> L(Phi,w+1) as y -> 0+, for all twelve rows.
   Theta(iy) = P(y) + (exponentially small), P a polynomial of degree w whose
   constant term is L(Phi,w+1); Richardson in y removes P's higher terms.   */
read("common.gp");
default(realprecision, 70);
Y0 = 0.004;
MM = 26000;
Th(cv,y) = my(u=exp(-2*Pi*y), M=min(MM,ceil(70*log(10)/(2*Pi*y))), s=0.); forstep(m=M,1,-1, s=(s+cv[m])*u); s;
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], w=r-1);
  my(cv=vector(MM,m,cm(nm,m)*1.0/m^(w+1)));
  my(f0=Th(cv,Y0), f1=Th(cv,Y0/2), f2=Th(cv,Y0/4), L=Ltarget(nm), ex);
  ex = if(w==1, 2*f1-f0, (8/3)*f2-2*f1+(1/3)*f0);
  print(nm,"\t Theta(0)=",ex,"\n\t L(Phi,",w+1,")=",L,"\n\t diff=",ex-L);
);
}
quit;
