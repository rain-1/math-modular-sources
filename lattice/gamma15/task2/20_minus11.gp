/* 20_minus11.gp -- conceptual origin of pi_D = -11*zeta(2)/5.
   OUTER orientation (1,psi): Phi_D = sum_n (sum_{d|n} w(d) d^2) q^n with
     w = re4 - 2*im4  (period 5, mean 0), so  D^{-2}Phi_D = sum_d w(d) Li_2(q^d).
   INNER orientation (psi,1): Phi = sum_n (sum_{d|n} v(n/d) d^2) q^n, so
     D^{-2}Phi = sum_k (v(k)/k^2) * q^k/(1-q^k).                              */
default(realprecision, 60);
re4(n) = my(r=n%5); if(r==1, 1, if(r==4, -1, 0));
im4(n) = my(r=n%5); if(r==2, 1, if(r==3, -1, 0));
s5=sqrt(5); ph5=(11+5*s5)/2; phm5=(5*s5-11)/2;
w(d) = re4(d) - 2*im4(d);
print("w(1..10) = ", vector(10,i,w(i)), "   (period 5, sum over a period = ", sum(a=1,5,w(a)), ")");
/* Abel regularisation of a P-periodic mean-zero sequence: lim sum a(d)x^d = (sum j a(j))/(-P) */
reg(f, P) = sum(j=1,P, j*f(j))/(-P);
print("");
print("=== OUTER: Li_2 structure.  cusp value = zeta(2) * reg(w * lambda) ===");
print("  cusp 0   (q -> +1): lambda_d = 1 for all d.  reg = ", reg(w, 5), "   -> pi = zeta(2)/5   [MATCHES]");
g(d) = if(d%2==0, 1, -1/2);   /* Li_2((-1)^d) / zeta(2) */
wg(d) = w(d)*g(d);
print("  wg(1..10) = ", vector(10,i,wg(i)), "   sum over period 10 = ", sum(a=1,10,wg(a)));
print("  cusp 1/2 (q -> -1): lambda_d = 1 (d even), -1/2 (d odd).  reg = ", reg(wg, 10), "  -> pi = -11*zeta(2)/5   [MATCHES]");
print("  sum_{j=1}^{10} j*wg(j) = ", sum(j=1,10,j*wg(j)), " ,  so reg = 22/(-10) = -11/5.");
print("");
print("=== INNER: Lambert structure.  1/(e^t-1) = 1/t - 1/2 + t/12 - ... ===");
v(k)  = 2*re4(k) - 2*ph5*im4(k);     /* Phi_new  */
vp(k) = 2*re4(k) + 2*phm5*im4(k);    /* Phi'_new = sigma(Phi_new) */
Lc(s) = 5^(-s)*(zetahurwitz(s,1/5)+I*zetahurwitz(s,2/5)-I*zetahurwitz(s,3/5)-zetahurwitz(s,4/5));
L2=Lc(2); L3=Lc(3); R2=real(L2); I2=imag(L2); R3=real(L3); I3=imag(L3);
print("  Re L(3,psi4) = ", R3, "    Im L(3,psi4) = ", I3);
print("  Re L(3,psi4) / Im L(3,psi4) = ", R3/I3, "     phi^5 = ", ph5);
print("  => sum_k v(k)/k^3 = 2(ReL3 - phi^5 ImL3) = ", 2*(R3-ph5*I3), "  (the cusp-0 convergence condition)");
print("  => sum_{k even} v'(k)/k^3 = (1/4)(-2 ImL3 + 2 phi^-5 ReL3) = ", (1/4)*(-2*I3+2*phm5*R3), "  (cusp-1/2 condition)");
print("     [both vanish <=> ReL(3,psi4) = phi^5 ImL(3,psi4): the SAME identity]");
print("  cusp 0,   Phi_new : -(1/2) sum_k v(k)/k^2  = ", -(R2 - ph5*I2), "     xi  = ", ph5*I2-R2);
print("  cusp 1/2, Phi'_new: -(1/2) sum_k v'(k)/k^2 = ", -(R2 + phm5*I2), "     xi' = ", -phm5*I2-R2);
print("  the weight -1/2 is the SAME for every k at BOTH cusps  =>  period = -(1/2) L(2,v),");
print("  which is K-linear in the source's K-coefficients  =>  exact sigma-equivariance.");
print("");
print("  (contrast: for Phi_new at cusp 1/2 the condition would be ReL3 = -phi^-5 ImL3 : ", 2*(R3+ph5*R3*0-ph5*I3)*0 + (-2*I3-2*ph5*R3)/1, " != 0, so Phi_new diverges there -- the log^2.)");
print("  11 = phi^5 - phi^-5 = L_5 (Lucas) = -(t1+t2) = Tr_{K/Q}(phi^5).");
quit;
