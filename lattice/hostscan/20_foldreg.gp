/* fold-regular subspace of M_3^Eis(Gamma_1(5)) on Zagier D's Apery-perfect host,
   and the periods it carries.  Basis Phi1..Phi4 (outer/inner, conjugate pairs). */
default(parisizemax, 10000000000);
default(realprecision, 140);
NQ = 110; NA = 84;
ps4(n) = if(n%5==0, 0, if(n%5==1, 1, if(n%5==2, I, if(n%5==3, -I, -1))));
ps4b(n) = conj(ps4(n));
one(n) = 1;
leg5(n) = kronecker(n,5);
xser(nq) = q*prod(n=1, nq-1, (1 - q^n + O(q^nq))^(5*leg5(n)));
Ap(n) = sum(k=0,n, binomial(n,k)^2*binomial(n+k,k));
eis(ch1, ch2, nq) = sum(n=1, nq-1, sumdiv(n, d, ch1(n/d)*ch2(d)*d^2)*q^n) + O(q^nq);
Dinv2(s) = my(n=serprec(s,q)); sum(k=1,n-1, polcoeff(s,k)/k^2*q^k) + O(q^n);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
xs = xser(NQ);
Fs = sum(n=0, NA, Ap(n)*xs^n) + O(q^NQ);
A = peel2(Fs, xs, NA, NQ);
PP = [eis(one,ps4,NQ), eis(one,ps4b,NQ), eis(ps4,one,NQ), eis(ps4b,one,NQ)];
/* constant terms of E_3^{psi,phi} at infinity: nonzero only when psi = 1 */
Bgen(ch) = -sum(a=1,5, ch(a)*(a^3/5 - 3*a^2/2 + 5*a/2*1 - 0))/1;
UU = vector(4, j, my(Th=Dinv2(PP[j])); my(B=peel2(Fs*Th, xs, NA, NQ)); vector(NA+1, i, B[i]*1.0/A[i]));
print("u_n for the four directions at n = NA, NA-1, NA-2:");
for(j=1,4, print("  Phi",j,": ", UU[j][NA+1], "  ", UU[j][NA], "  ", UU[j][NA-1]));
/* differences: the non-fold-regular part */
D1 = vector(4, j, UU[j][NA+1]-UU[j][NA]);
D2 = vector(4, j, UU[j][NA]-UU[j][NA-1]);
D3 = vector(4, j, UU[j][NA-1]-UU[j][NA-2]);
print("differences d_n: ", D1);
print("               : ", D2);
print("ratios d_n/d_{n-1}: ", vector(4,j, D1[j]/D2[j]));
M = matrix(3,4,i,j, if(i==1,D1[j], if(i==2,D2[j], D3[j])));
K = matker(M);
print("kernel dim of the difference matrix (numerical): ", matsize(K));
/* explicit rational candidates */
tst(cv, nm) = my(Ph = sum(j=1,4, cv[j]*PP[j])); my(Th=Dinv2(Ph)); my(B=peel2(Fs*Th, xs, NA, NQ)); print(nm,":  xi_n = ", B[NA+1]*1.0/A[NA+1], "   xi_{n-1} = ", B[NA]*1.0/A[NA]);
tst([1,1,0,0], "Phi1+Phi2 (outer real)");
tst([0,0,1,1], "Phi3+Phi4 (inner real)");
tst([1,1,-1,-1], "outer - inner");
tst([1,1,1,1], "outer + inner");
tst([-2,-2,1,1], "-2*outer + inner");
tst([2,2,1,1], "2*outer + inner");
Lchi(ch, s) = 5^(-s)*sum(a=1,4, ch(a)*zetahurwitz(s, a/5));
print("zeta(2) = ", zeta(2));
print("Re L(2,psi4) = ", real(Lchi(ps4,2)), "   Im = ", imag(Lchi(ps4,2)));
quit;
