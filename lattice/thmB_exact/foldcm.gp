/* foldcm.gp -- the interior (Fricke) fold of the four real third-order rows:
   tau_* = i/sqrt(N) is the fold, Phi(tau_*)=0, L(Phi,2)=0, Fricke sign -1,
   and the fold-connection value equals L(Phi,3).                          */
read("common.gp");
default(realprecision, 70);
MM = 20000;
NLEV = Map(); mapput(NLEV,"alpha",12); mapput(NLEV,"gamma",6); mapput(NLEV,"eps",8);
mapput(NLEV,"zeta",9); mapput(NLEV,"delta",12); mapput(NLEV,"eta",20);
NT = 700;
ev(cv,k,q0) = my(M=min(MM,ceil(80*log(10)/(-log(abs(q0))))), s=0.); forstep(m=M,1,-1, s=s*q0 + cv[m]*m^k); s*q0;
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], a=R[3], b=R[4], c=R[5], w=r-1);
  if(r!=3, next);
  my(N=mapget(NLEV,nm));
  my(cv=vector(MM,m,cm(nm,m)*1.0));
  my(cvT=vector(MM,m,cv[m]/m^(w+1)));
  my(tau=I/sqrt(N), q0=exp(2*Pi*I*tau));
  my(BB=build(r,a,b,c,NT), tq=truncate(BB[1]), Fq=truncate(BB[2]));
  my(Fv=subst(Fq,'t,q0), Fd=subst(deriv(Fq,'t),'t,q0));
  my(td=subst(deriv(tq,'t),'t,q0), tv=subst(tq,'t,q0));
  my(Th=ev(cvT,0,q0)/q0*q0, Thd=0.);
  \\ Theta(q)=sum c(m) m^-(w+1) q^m ; Theta'(q)= sum m c(m) m^-(w+1) q^(m-1)
  Th  = ev(cvT,0,q0);
  Thd = ev(cvT,1,q0)/q0;
  my(xi = Th + Fv*Thd/Fd, L3=Ltarget(nm));
  print("=== ",nm,"  N=",N,"  tau_*=i/sqrt(",N,")  q_c=",q0);
  print("   t(q_c)      = ", tv, "     (root of P: ",polroots(1-2*a*x+c*x^2)~,")");
  print("   t'(q_c)     = ", td);
  print("   Phi(q_c)    = ", ev(cv,0,q0));
  print("   F(q_c)      = ", Fv, "   F'(q_c)=",Fd);
  print("   xi (fold)   = ", xi);
  print("   L(Phi,3)    = ", L3);
  print("   xi - L      = ", xi-L3);
  print("   L(Phi,2)    = ", if(nm=="zeta", lfun(-3,2)*lfun(-3,-1), Pmellin(nm,2)*zeta(2)*zeta(-1)));
  print("   P(w+2)=P(4) = ", Pmellin(nm,4));
);
}
quit;
