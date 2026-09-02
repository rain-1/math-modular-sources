/* Dump A_n, B_{D,n} and B_{new,n} = (u+v sqrt5) for the far-cusp period computation. */
default(parisizemax, 12000000000);
NQ = 260; NA = 220;
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
xs = xser(NQ); Fs = sum(n=0, NA, Ap(n)*xs^n) + O(q^NQ); A = peel2(Fs, xs, NA, NQ);
P1=eis(one,ps4,NQ); P2=eis(one,ps4b,NQ); P3=eis(ps4,one,NQ); P4=eis(ps4b,one,NQ);
R1 = realser(P1+P2, NQ); R2 = realser(I*(P1-P2), NQ);
R3 = realser(P3+P4, NQ); R4 = realser(I*(P3-P4), NQ);
PhiD = R1/2 + R2;
BD = peel2(Fs*Dinv2(PhiD), xs, NA, NQ);
/* Phi_new = R3 + phi^5 R4 with phi^5 = (11+5 sqrt5)/2 : rational and sqrt5 parts */
B3 = peel2(Fs*Dinv2(R3), xs, NA, NQ);
B4 = peel2(Fs*Dinv2(R4), xs, NA, NQ);
/* B_new = B3 + (11/2) B4 + (5/2) sqrt5 B4 ; B'_new = B3 + (11/2) B4 - (5/2) sqrt5 B4 */
write("farcusp_data.txt", "# n  A_n  BD_n  U_n  V_n   with B_new = U + V sqrt5, B'_new = U - V sqrt5");
for(n=0, NA, write("farcusp_data.txt", n, " ", A[n+1], " ", BD[n+1], " ", B3[n+1]+11/2*B4[n+1], " ", 5/2*B4[n+1]));
print("wrote farcusp_data.txt, ", NA+1, " rows");
print("A     = ", vector(6,i,A[i]));
print("BD    = ", vector(6,i,BD[i]));
print("U     = ", vector(6,i,B3[i]+11/2*B4[i]));
print("V     = ", vector(6,i,5/2*B4[i]));
quit;
