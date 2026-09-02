/* Independent re-verification of the new period xi on Zagier's row D / Gamma_1(5).
   xi = lim B_n/A_n for the fold-regular inner Eisenstein direction
        Phi_new = R3 + phi^5 R4,  R3 = E3(psi4,1)+E3(psi4bar,1), R4 = i(E3(psi4,1)-E3(psi4bar,1)).
   Also re-derives Phi_D = R1/2 + R2 with Apery limit zeta(2)/5. */
default(parisizemax, 12000000000);
default(realprecision, 160);
NQ = 150; NA = 118;
ps4(n) = if(n%5==0, 0, if(n%5==1, 1, if(n%5==2, I, if(n%5==3, -I, -1))));
ps4b(n) = conj(ps4(n));
one(n) = 1;
leg5(n) = kronecker(n,5);
xser(nq) = q*prod(n=1, nq-1, (1 - q^n + O(q^nq))^(5*leg5(n)));
Ap(n) = sum(k=0,n, binomial(n,k)^2*binomial(n+k,k));
eis(c1, c2, nq) = sum(n=1, nq-1, sumdiv(n, d, c1(n/d)*c2(d)*d^2)*q^n) + O(q^nq);
Dinv2(s) = my(n=serprec(s,q)); sum(k=1,n-1, polcoeff(s,k)/k^2*q^k) + O(q^n);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
realser(s, nq) = sum(k=0,nq-1, real(polcoeff(s,k))*q^k) + O(q^nq);
xs = xser(NQ);
Fs = sum(n=0, NA, Ap(n)*xs^n) + O(q^NQ);
A = peel2(Fs, xs, NA, NQ);
print("A_n (row D) = ", vector(8,i,A[i]));
P1 = eis(one,ps4,NQ); P2 = eis(one,ps4b,NQ); P3 = eis(ps4,one,NQ); P4 = eis(ps4b,one,NQ);
R1 = realser(P1+P2, NQ); R2 = realser(I*(P1-P2), NQ);
R3 = realser(P3+P4, NQ); R4 = realser(I*(P3-P4), NQ);
ph5 = (11+5*sqrt(5))/2;
PhiD = R1/2 + R2;
PhiN = R3 + ph5*R4;
print("Phi_D  q-coeffs = ", vector(10,i,polcoeff(PhiD,i-1)));
print("Phi_new q-coeffs (numeric) = ", vector(6,i,polcoeff(PhiN,i-1)));
BD = peel2(Fs*Dinv2(PhiD), xs, NA, NQ);
BN = peel2(Fs*Dinv2(PhiN), xs, NA, NQ);
uD1 = BD[NA+1]/A[NA+1]; uD2 = BD[NA]/A[NA];
uN1 = BN[NA+1]/A[NA+1]; uN2 = BN[NA]/A[NA];
print("zeta(2)/5              = ", zeta(2)/5);
print("B_D,n / A_n  (n=", NA, ")  = ", uD1);
print("  agreement digits n vs n-1 : ", floor(-log(abs(uD1-uD2)+1e-300)/log(10)));
print("  |B_D,n/A_n - zeta(2)/5|   = ", abs(uD1 - zeta(2)/5));
print("xi (n=", NA, ")   = ", uN1);
print("xi (n=", NA-1, ") = ", uN2);
print("  agreement digits n vs n-1 : ", floor(-log(abs(uN1-uN2)+1e-300)/log(10)));
Lchi(ch, s) = 5^(-s)*sum(a=1,4, ch(a)*zetahurwitz(s, a/5));
LR = real(Lchi(ps4,2)); LI = imag(Lchi(ps4,2));
print("Re L(2,psi4) = ", LR);
print("Im L(2,psi4) = ", LI);
print("phi^5 Im L - Re L = ", ph5*LI - LR);
print("  |xi - (phi^5 Im L - Re L)| = ", abs(uN1 - (ph5*LI - LR)));
print("lindep [xi,1,zeta(2),ReL,ImL,sqrt5,sqrt5*ReL,sqrt5*ImL] = ", lindep([uN1,1,zeta(2),LR,LI,sqrt(5),sqrt(5)*LR,sqrt(5)*LI],40));
quit;
