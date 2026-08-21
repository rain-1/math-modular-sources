/* cvals.gp -- the complex-fold connection constants of delta and eta, to high
   precision, from the fold-connection formula at the exact fold point.    */
read("common.gp");
default(realprecision, 90);
MM = 30000;
NT = 900;
ev(cv,k,q0,pr) = my(M=min(MM,ceil(pr*log(10)/(-log(abs(q0))))), s=0.); forstep(m=M,1,-1, s=s*q0 + cv[m]*m^k); s*q0;
{
DATA = [["delta",3,7,3,81, 0.17+0.24*I], ["eta",3,11,5,125, 0.2+0.1*I], ["eta",3,11,5,125, 0.1+0.2*I]];
for(i=1,#DATA,
  my(D=DATA[i], nm=D[1], r=D[2], a=D[3], b=D[4], c=D[5], tau0=D[6], w=r-1);
  my(BB=build(r,a,b,c,NT), tq=truncate(BB[1]), Fq=truncate(BB[2]), dt=deriv(tq,'t), d2t=deriv(dt,'t));
  my(q0=exp(2*Pi*I*tau0));
  for(it=1,200, my(v=subst(dt,'t,q0), dv=subst(d2t,'t,q0)); q0 = q0 - v/dv);
  my(tau=log(q0)/(2*Pi*I));
  my(cv=vector(MM,m,cm(nm,m)*1.0), cvT=vector(MM,m,cm(nm,m)*1.0/m^(w+1)));
  my(Fv=subst(Fq,'t,q0), Fd=subst(deriv(Fq,'t),'t,q0));
  my(Th=ev(cvT,0,q0,95), Thd=ev(cvT,1,q0,95)/q0);
  my(xi=Th+Fv*Thd/Fd, L=Ltarget(nm));
  print("=== ",nm,"  tau_* = ",tau);
  print("   |tau|^2*N? : 12|tau|^2=",12*abs(tau)^2,"  20|tau|^2=",20*abs(tau)^2);
  print("   q_c   = ",q0,"   |q_c|=",abs(q0));
  print("   t(q_c)= ",subst(tq,'t,q0));
  print("   t'(q_c)=",subst(dt,'t,q0));
  print("   Phi(q_c)=",ev(cv,0,q0,95));
  print("   F(q_c)= ",Fv);
  print("   xi    = ",xi);
  print("   L(Phi,",w+1,") = ",L);
  print("   xi/L  = ",xi/L);
  print("   xi - L= ",xi-L);
);
}
quit;
