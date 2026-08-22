/* a6_design.gp -- the design-rule / master-formula numbers for the L(3,chi5) class at p=5. */
default(realprecision, 30);
L1e = sqrt(125.); L1f = 20*sqrt(5.);      /* |lambda| for eta and for AESZ 184 */
k = 3; s5 = 3; l5 = log(5.);
print("eta : log Lambda = ", log(L1e), "   AESZ184: log Lambda = ", log(L1f));
print("3*log5 = ", 3*l5);
budget(L,kk,sig) = -log(L) - kk + sig*l5;   /* log(1/|lambda_2|) - k + sum sigma_p log p, lambda_2=lambda_1 */
print("\nbudget(eta)    = ", budget(L1e,k,s5));
print("budget(AESZ184)= ", budget(L1f,k,s5));
print("(Hadamard twist costs exactly log 4 = ", log(4.), ")");
/* master formula, no decay: lambda = Lambda */
F(Ldec,Leng,kk,sig,r) = 0.5*( kk*max(r,1) + r*log(Leng) - ( -log(Ldec) + min(sig*r,sig)*l5 ) );
print("\n--- master formula, decayer=184 engine=eta, r=alpha/gamma ---");
forstep(r=0.5,3,0.25, my(f=F(L1f,L1e,k,s5,r)); print("  r=",r,"  F=",f,"  H=F+0 => delta=0"));
print("--- decayer=eta engine=184 ---");
forstep(r=0.5,3,0.25, my(f=F(L1e,L1f,k,s5,r)); print("  r=",r,"  F=",f));
print("\n--- threshold for a hypothetical decayer against engine eta ---");
print("  engine eta: k_eng=3, sigma_5=3, log rho_2^eng = ", log(L1e));
print("  at balance r = sigma_dec/3, F<0  <=>  log(1/lambda) > max(sigma_dec,k_d) + nu + kappa_5*log5 - (log(sqrt125)/3)*sigma_dec");
print("  coefficient (log sqrt125)/3 = ", log(L1e)/3, " ;  1 - that = ", 1-log(L1e)/3);
print("  so with k_d <= sigma_dec:  log(1/lambda) > ", 1-log(L1e)/3, "*sigma_dec + nu + ", l5, "*kappa_5");
print("  and sigma_dec = w + 2 kappa_5 =>  log(1/lambda) > ", 1-log(L1e)/3, "*w + ", 2*(1-log(L1e)/3)+l5, "*kappa_5 + nu");
quit
