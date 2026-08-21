/* 03_fricke_fold.gp -- the six third-order rows.
   (i) Fricke data: Phi|_4 W_N = e Phi, F|_2 W_N = e F, t o W_N = t;
   (ii) the fold: tau_* = i/sqrt(N), t'(q_c)=0, t(q_c)=t_c, Phi(q_c)=0;
   (iii) the fold-connection value xi = Theta + F Theta'/F' vs L(Phi,3);
   (iv) the middle critical value L(Phi,2).                                */
read("common.gp");
default(realprecision, 70);
MM = 20000; NT = 700;
NLEV = Map(); mapput(NLEV,"alpha",12); mapput(NLEV,"gamma",6); mapput(NLEV,"eps",8);
mapput(NLEV,"zeta",9); mapput(NLEV,"delta",12); mapput(NLEV,"eta",20);
ev(cv,k,q0) = my(M=min(MM,ceil(80*log(10)/(-log(abs(q0))))), s=0.); forstep(m=M,1,-1, s=s*q0 + cv[m]*m^k); s*q0;
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], a=R[3], b=R[4], c=R[5], w=r-1);
  if(r!=3, next);
  my(N=mapget(NLEV,nm));
  my(cv=vector(MM,m,cm(nm,m)*1.0), cvT=vector(MM,m,cm(nm,m)*1.0/m^(w+1)));
  my(tau=I/sqrt(N), q0=exp(2*Pi*I*tau));
  my(BB=build(r,a,b,c,NT), tq=truncate(BB[1]), Fq=truncate(BB[2]));
  my(Fv=subst(Fq,'t,q0), Fd=subst(deriv(Fq,'t),'t,q0));
  my(Th=ev(cvT,0,q0), Thd=ev(cvT,1,q0)/q0);
  my(xi=Th+Fv*Thd/Fd, L3=Ltarget(nm));
  print("=== ",nm,"  N=",N);
  my(tt=0.19+0.53*I, ttp=-1/(N*tt), Q0=exp(2*Pi*I*tt), Q1=exp(2*Pi*I*ttp), j=sqrt(N)*tt);
  print("   Phi|_4 W_N / Phi = ", ev(cv,0,Q1)/(j^4*ev(cv,0,Q0)));
  print("   F|_2 W_N   / F   = ", subst(Fq,'t,Q1)/(j^2*subst(Fq,'t,Q0)));
  print("   t(W_N tau) - t(tau) = ", subst(tq,'t,Q1)-subst(tq,'t,Q0));
  print("   t(q_c)-t_c  = ", subst(tq,'t,q0)-vecsort(polroots(1-2*a*x+c*x^2),abs)[1]);
  print("   t'(q_c)     = ", subst(deriv(tq,'t),'t,q0));
  print("   Phi(q_c)    = ", ev(cv,0,q0));
  print("   xi(fold)-L(Phi,3) = ", xi-L3, "     [L(Phi,3) = ",L3,"]");
  print("   L(Phi,2)    = ", if(nm=="D",0,LPhi(nm,2)));
  print("   endpoint delta(phi)P(4) = ", endpoint(nm,r));
);
}
quit;
