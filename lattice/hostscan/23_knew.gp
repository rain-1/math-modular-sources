/* denominator exponent of the new linear form on Gamma_1(5) */
default(parisizemax, 10000000000);
default(realprecision, 60);
NQ = 100; NA = 78;
ps4(n) = if(n%5==0, 0, if(n%5==1, 1, if(n%5==2, I, if(n%5==3, -I, -1))));
ps4b(n) = conj(ps4(n));
one(n) = 1;
leg5(n) = kronecker(n,5);
xser(nq) = q*prod(n=1, nq-1, (1 - q^n + O(q^nq))^(5*leg5(n)));
Ap(n) = sum(k=0,n, binomial(n,k)^2*binomial(n+k,k));
eisI(ch1, ch2, nq) = sum(n=1, nq-1, sumdiv(n, d, ch1(n/d)*ch2(d)*d^2)*q^n) + O(q^nq);
Dinv2(s) = my(n=serprec(s,q)); sum(k=1,n-1, polcoeff(s,k)/k^2*q^k) + O(q^n);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
xs = xser(NQ); Fs = sum(n=0, NA, Ap(n)*xs^n) + O(q^NQ);
P3=eisI(ps4,one,NQ); P4=eisI(ps4b,one,NQ);
R3 = P3+P4; R4 = I*(P3-P4);
/* work over Z[w], w = sqrt5:  phi^5 = (11+5w)/2 ; write cand = R3 + phi5*R4 with w symbolic */
w = Mod(x, x^2-5);
ph5 = (11+5*w)/2;
cand = R3 + ph5*R4;
print("cand coefficients (in Z[(1+sqrt5)/2]) c(1..10) = ", vector(10,i,polcoeff(cand,i-1)));
Th = Dinv2(cand);
yB = Fs*Th;
B = peel2(yB, xs, NA, NQ);
print("B_1..B_6 = ", vector(6,i,B[i+1]));
/* denominator exponent: d_n^k * B_n integral in Z[w]? */
dn = 1; km = 0;
for(n=1, NA, dn = lcm(dn, n); my(bn = B[n+1]); if(bn != 0, my(v = liftall(bn)); my(de = lcm(denominator(polcoeff(v,0)), denominator(polcoeff(v,1)))); my(kj=0, t=de); while(t>1 && kj<12, kj++; t = t/gcd(t,dn)); if(kj>km, km=kj; print("   n=",n," k=",kj," den=",de))));
print("denominator exponent k = ", km, "   (Zagier D itself has k = 2)");
quit;
