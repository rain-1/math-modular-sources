read("common.gp");
default(realprecision, 60);
MM = 4000;
CC = Map();
getc(nm) = { my(v=vector(MM)); for(m=1,MM, v[m]=cm(nm,m)); v; }
evalPhi(cv, q0, nd) = { my(M=ceil(nd*log(10)/(-log(abs(q0)))), s=0.); M=min(M,MM);
   forstep(m=M,1,-1, s = (s + cv[m])*q0); s; }
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], a=R[3], b=R[4], c=R[5]);
  my(P=if(r==2,1-a*x+c*x^2,1-2*a*x+c*x^2), rts=polroots(P));
  if(abs(imag(rts[1]))<1e-20, next);
  my(cv=getc(nm), best=List());
  print("=== ",nm," (complex fold) t_c=",rts~);
  for(ix=-50,50, for(iy=4,50,
    my(tau=ix/100. + I*iy/100., q0=exp(2*Pi*I*tau), v=evalPhi(cv,q0,70));
    if(abs(v) < 0.02*abs(q0), listput(best,tau))));
  print("  candidate tau (coarse |Phi| small): ", Vec(best));
);
}
quit;
