/* 04_complex_folds.gp -- the three complex-fold rows.

   delta, eta: the fold is an interior CM point on the Fricke circle
   |tau| = 1/sqrt(N); the connection value is Theta + F Theta'/F' there.
   B: t_c sits at the cusp 1/6 of Gamma_0(36); with A = [[1,0],[6,1]] one has
   Theta(A sigma)(6 sigma + 1) = c_1 sigma + c_0 + O(exp(-2 pi Im sigma)),
   and the connection value is c_1/6.                                      */
read("common.gp");
default(realprecision, 90);
MM = 900000; NT = 900;
ev(cv,k,q0,pr) = my(M=min(MM,ceil(pr*log(10)/(-log(abs(q0))))), s=0.); forstep(m=M,1,-1, s=s*q0 + cv[m]*m^k); s*q0;

print("======== interior complex folds: delta, eta ========");
{
DATA = [["delta",3,7,3,81, 0.17+0.24*I],["eta",3,11,5,125, 0.1+0.2*I],["eta",3,11,5,125, 0.2+0.1*I]];
for(i=1,#DATA,
  my(D=DATA[i], nm=D[1], r=D[2], a=D[3], b=D[4], c=D[5], tau0=D[6], w=r-1);
  my(BB=build(r,a,b,c,NT), tq=truncate(BB[1]), Fq=truncate(BB[2]),
     dt=deriv(tq,'t), d2t=deriv(dt,'t));
  my(q0=exp(2*Pi*I*tau0));
  for(it=1,200, q0 = q0 - subst(dt,'t,q0)/subst(d2t,'t,q0));
  my(tau=log(q0)/(2*Pi*I), N=if(nm=="delta",12,20));
  my(cvT=vector(MM,m,cm(nm,m)*1.0/m^(w+1)));
  my(Fv=subst(Fq,'t,q0), Fd=subst(deriv(Fq,'t),'t,q0));
  my(xi = ev(cvT,0,q0,95) + Fv*(ev(cvT,1,q0,95)/q0)/Fd, L=Ltarget(nm));
  print("--- ",nm,"  tau_* = ",tau,"    N|tau_*|^2 = ",N*abs(tau)^2);
  print("    t(q_c)  = ",subst(tq,'t,q0),"   t'(q_c)=",subst(dt,'t,q0));
  print("    xi      = ",xi);
  print("    Re xi - L(Phi,3) = ",real(xi)-L);
  if(nm=="delta",
     print("    Im xi - Pi^3/81  = ",imag(xi)-Pi^3/81);
     print("    Im xi / L(Phi,2) = ",imag(xi)/LPhi(nm,2),"     2*Pi/3 = ",2*Pi/3),
     print("    Im xi / (Pi*L(chi5,2)) = ",imag(xi)/(Pi*lfun(5,2))));
);
}

print("");
print("======== row B: the cusp 1/6 of Gamma_0(36) ========");
{
default(realprecision,60);
my(cv = vector(MM, m, cm("B",m)*1.0/m^2));
G(Y) = my(sig=I*Y, tau=sig/(6*sig+1), q0=exp(2*Pi*I*tau),
          M=min(MM,ceil(55*log(10)/(-log(abs(q0))))), s=0.);
        forstep(m=M,1,-1, s=(s+cv[m])*q0); s*(6*sig+1);
my(Y1=18., Y2=23., g1=G(Y1), g2=G(Y2));
my(c1=(g1-g2)/(I*Y1-I*Y2), xiB=c1/6);
print("    xi_B = ",xiB);
print("    Re xi_B - L(Phi_B,2)        = ",real(xiB)-Ltarget("B"));
print("    Im xi_B - 2*Pi^2/(27*sqrt3) = ",imag(xiB)-2*Pi^2/(27*sqrt(3)));
print("    lindep(Im xi_B, Pi^2/sqrt3) = ",lindep([imag(xiB),Pi^2/sqrt(3)]));
}
quit;
