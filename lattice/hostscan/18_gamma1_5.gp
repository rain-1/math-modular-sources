/* 18_gamma1_5.gp -- Zagier D's Apery-perfect host Gamma_1(5) (lam2^norm = 1, k = 2):
   which periods does the weight-3 Eisenstein source space carry?
   psi4 = odd quartic character mod 5 with psi4(2) = i.  Directions:
     Phi1 = E_3^{1,psi4}, Phi2 = E_3^{1,psi4bar}   (outer)  -> period zeta(2)
     Phi3 = E_3^{psi4,1}, Phi4 = E_3^{psi4bar,1}   (inner)  -> period L(2,psi4), L(2,psi4bar)
   Row: x = q prod (1-q^n)^{5 legendre(n,5)},  F = sum A_n x^n with A_n Apery's zeta(2) numbers. */
default(parisizemax, 8000000000);
default(realprecision, 90);
NQ = 200;
ps4(n) = if(n%5==0, 0, if(n%5==1, 1, if(n%5==2, I, if(n%5==3, -I, -1))));
ps4b(n) = conj(ps4(n));
leg5(n) = kronecker(n,5);
xser(nq) = q*prod(n=1, nq-1, (1 - q^n + O(q^nq))^(5*leg5(n)));
/* A_n = sum_k binom(n,k)^2 binom(n+k,k)  (Apery zeta(2)) */
Ap(n) = sum(k=0,n, binomial(n,k)^2*binomial(n+k,k));
eis(ch1, ch2, nq) = sum(n=1, nq-1, sumdiv(n, d, ch1(n/d)*ch2(d)*d^2)*q^n) + O(q^nq);
one(n) = 1;
Dinv2(s) = my(n=serprec(s,q)); sum(k=1,n-1, polcoeff(s,k)/k^2*q^k) + O(q^n);
xs = xser(NQ);
A = vector(60, i, Ap(i-1));
/* build F = sum A_n x^n */
Fs = sum(n=0, 55, A[n+1]*xs^n) + O(q^56);
print("F = ", vector(12,i,polcoeff(Fs,i-1)));
print("x = ", vector(12,i,polcoeff(xs,i-1)));
/* recover A_n from F and x by peeling, as a check */
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
print("peel check: ", vector(8,i,peel2(Fs,xs,20,56)[i]));
quit;
