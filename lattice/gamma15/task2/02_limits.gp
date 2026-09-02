/* 02_limits.gp -- Apery limits, closed forms, lindep certificates, decay rates. */
default(parisizemax, 20000000000);
default(realprecision, 700);
read("build.txt");
s5 = sqrt(5); ph5 = (11+5*s5)/2; phm5 = (5*s5-11)/2;   /* phi^5, phi^{-5} */
print("phi^5   = ", ph5);
print("phi^-5  = ", phm5);
print("phi^5*phi^-5 = ", ph5*phm5);
t1 = phm5; t2 = -ph5;
print("check 1-11t-t^2 at t1,t2: ", 1-11*t1-t1^2, "  ", 1-11*t2-t2^2);
Lchi(s) = 5^(-s)*(zetahurwitz(s,1/5) + I*zetahurwitz(s,2/5) - I*zetahurwitz(s,3/5) - zetahurwitz(s,4/5));
LL = Lchi(2); LR = real(LL); LI = imag(LL);
print("Re L(2,psi4) = ", LR);
print("Im L(2,psi4) = ", LI);
z2 = zeta(2);
print("zeta(2)/5    = ", z2/5);
xipred  = ph5*LI - LR;
xippred = -phm5*LI - LR;
print("xi  predicted = phi^5 ImL - ReL   = ", xipred);
print("xi' predicted = -phi^-5 ImL - ReL = ", xippred);
/* companions */
Bnew  = vector(NA+1, i, B3[i] + ph5*B4[i]);
Bnewp = vector(NA+1, i, B3[i] - phm5*B4[i]);
rat(n) = BD[n+1]/A[n+1];
print("");
print("=== Apery limits, n = NA and NA-1 ===");
for(j=0,2, my(n=NA-j); print("  BD_n/A_n   (n=",n,") = ", BD[n+1]*1.0/A[n+1]));
for(j=0,2, my(n=NA-j); print("  Bnew_n/A_n (n=",n,") = ", Bnew[n+1]/A[n+1]));
for(j=0,2, my(n=NA-j); print("  Bnew'_n/A_n(n=",n,") = ", Bnewp[n+1]/A[n+1]));
dig(a,b) = if(a==b, 9999, -log(abs(a-b))/log(10));
print("");
print("digits BD_n/A_n vs zeta(2)/5 : ", dig(BD[NA+1]*1.0/A[NA+1], z2/5));
print("digits Bnew_n/A_n vs xipred  : ", dig(Bnew[NA+1]/A[NA+1], xipred));
print("digits Bnew'_n/A_n vs xippred: ", dig(Bnewp[NA+1]/A[NA+1], xippred));
print("digits Bnew'_n/A_n vs Bnew'_{n-1}/A_{n-1}: ", dig(Bnewp[NA+1]/A[NA+1], Bnewp[NA]/A[NA]));
print("digits BD_n/A_n vs BD_{n-1}/A_{n-1}: ", dig(BD[NA+1]*1.0/A[NA+1], BD[NA]*1.0/A[NA]));
print("digits Bnew_n/A_n vs Bnew_{n-1}/A_{n-1}: ", dig(Bnew[NA+1]/A[NA+1], Bnew[NA]/A[NA]));
quit;
