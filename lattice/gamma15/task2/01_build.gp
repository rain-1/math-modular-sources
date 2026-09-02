/* 01_build.gp -- exact construction of the Gamma_1(5) Zagier-D host, the two
   weight-3 Eisenstein sources and their Apery companions, all over Q.
   Writes A[n], BD[n], B3[n], B4[n] (exact rationals) to build.txt.
   B_new = B3 + phi^5 B4 ; B_new' = sigma(B_new) = B3 - phi^{-5} B4.        */
default(parisizemax, 30000000000);
NA = 200;                      /* highest coefficient index                 */
NQ = NA + 3;                   /* q-adic precision                          */
leg5(n) = kronecker(n,5);
re4(n) = my(r=n%5); if(r==1, 1, if(r==4, -1, 0));
im4(n) = my(r=n%5); if(r==2, 1, if(r==3, -1, 0));
Ap(n) = sum(k=0, n, binomial(n,k)^2*binomial(n+k,k));
xser(nq) = q*prod(n=1, nq-1, (1 - q^n + O(q^nq))^(5*leg5(n)));
Dinv2(s) = my(m=serprec(s,q)); sum(k=1, m-1, polcoeff(s,k)/k^2*q^k) + O(q^m);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
gettime();
xs = xser(NQ);
print("x built, ", gettime(), " ms;  x = ", vector(9,i,polcoeff(xs,i-1)));
Fs = sum(n=0, NQ-1, Ap(n)*xs^n) + O(q^NQ);
print("F built, ", gettime(), " ms;  F = ", vector(9,i,polcoeff(Fs,i-1)));
A = peel2(Fs, xs, NA, NQ);
print("peel check A_0..A_7 = ", vector(8,i,A[i]), "  (", gettime(), " ms)");
/* the four real Eisenstein directions, rational q-expansions */
X = vector(NQ-1, n, sumdiv(n, d, re4(d)*d^2));
Y = vector(NQ-1, n, sumdiv(n, d, im4(d)*d^2));
U = vector(NQ-1, n, sumdiv(n, d, re4(n/d)*d^2));
V = vector(NQ-1, n, sumdiv(n, d, im4(n/d)*d^2));
mkser(cv) = sum(n=1, NQ-1, cv[n]*q^n) + O(q^NQ);
R1 = 2*mkser(X); R2 = -2*mkser(Y); R3 = 2*mkser(U); R4 = -2*mkser(V);
PhiD = R1/2 + R2;
print("Phi_D  = ", vector(9,i,polcoeff(PhiD,i-1)));
print("R1 = ", vector(9,i,polcoeff(R1,i-1)));
print("R2 = ", vector(9,i,polcoeff(R2,i-1)));
print("R3 = ", vector(9,i,polcoeff(R3,i-1)));
print("R4 = ", vector(9,i,polcoeff(R4,i-1)));
BD = peel2(Fs*Dinv2(PhiD), xs, NA, NQ);
print("BD done ", gettime(), " ms");
B3 = peel2(Fs*Dinv2(R3), xs, NA, NQ);
print("B3 done ", gettime(), " ms");
B4 = peel2(Fs*Dinv2(R4), xs, NA, NQ);
print("B4 done ", gettime(), " ms");
write("build.txt", "NA = ", NA);
write("build.txt", "A = ", A);
write("build.txt", "BD = ", BD);
write("build.txt", "B3 = ", B3);
write("build.txt", "B4 = ", B4);
print("BD_0..BD_5 = ", vector(6,i,BD[i]));
print("B3_0..B3_5 = ", vector(6,i,B3[i]));
print("B4_0..B4_5 = ", vector(6,i,B4[i]));
quit;
