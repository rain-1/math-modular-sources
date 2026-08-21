/* Bcusp.gp -- row B: t_c is attained at the cusp 1/6 of Gamma_0(36); the
   connection constant is lim Theta(tau) as tau -> 1/6 along A(i Y),
   A=[[1,0],[6,1]].                                                       */
read("common.gp");
default(realprecision, 60);
MM = 400000;
cv = vector(MM, m, cm("B",m)*1.0/m^2);
Th(q0,pr) = my(M=min(MM,ceil(pr*log(10)/(-log(abs(q0))))), s=0.); forstep(m=M,1,-1, s=(s+cv[m])*q0); s;
{
NT=700; BB=build(2,9,3,27,NT); tq=truncate(BB[1]);
print("t at 1/6 + i*eps:");
for(k=1,5, my(eps=0.02/2^k, tau=1/6.+I*eps, q0=exp(2*Pi*I*tau));
   print("  eps=",eps,"  t=",subst(tq,'t,q0)));
print("t_c = ", polroots(1-9*x+27*x^2)~);
print("");
print("Theta along A(iY), A=[[1,0],[6,1]]:");
vals=List();
for(k=0,6, my(Y=10.*2^k, sig=I*Y, tau=sig/(6*sig+1), q0=exp(2*Pi*I*tau), v=Th(q0,55));
   print("  Y=",Y,"  Im tau=",imag(tau),"  Theta=",v); listput(vals,[Y,v]));
V=Vec(vals);
print("");
print("Richardson (limit = (2*f(2Y)-f(Y)) style in 1/Y):");
for(k=1,#V-1, print("   ",2*V[k+1][2]-V[k][2]));
print("");
print("L(Phi_B,2) = ", Ltarget("B"));
}
quit;
