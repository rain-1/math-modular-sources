/* 15_checks.gp -- calibration checks. */
default(parisizemax, 20000000000);
NQ = 40;
leg5(n) = kronecker(n,5);
re4(n) = my(r=n%5); if(r==1, 1, if(r==4, -1, 0));
im4(n) = my(r=n%5); if(r==2, 1, if(r==3, -1, 0));
Ap(n) = sum(k=0, n, binomial(n,k)^2*binomial(n+k,k));
xser(nq) = q*prod(n=1, nq-1, (1 - q^n + O(q^nq))^(5*leg5(n)));
xs = xser(NQ); Fs = sum(n=0, NQ-1, Ap(n)*xs^n) + O(q^NQ);
mk(f) = sum(n=1, NQ-1, f(n)*q^n) + O(q^NQ);
R1 = 2*mk(n->sumdiv(n,d,re4(d)*d^2)); R2 = -2*mk(n->sumdiv(n,d,im4(d)*d^2));
R3 = 2*mk(n->sumdiv(n,d,re4(n/d)*d^2)); R4 = -2*mk(n->sumdiv(n,d,im4(n/d)*d^2));
th(s) = my(m=serprec(s,q)); sum(k=1,m-1, k*polcoeff(s,k)*q^k) + O(q^m);
PhiD = R1/2 + R2;
print("Phi_D  - F*theta_q(x)  = ", PhiD - Fs*th(xs));
print("Phi_D  q-exp = ", vector(10,i,polcoeff(PhiD,i-1)));
/* Phi_new' = (1-I*phi^-5)P3 + (1+I*phi^-5)P4 ? */
P3 = mk(n->sumdiv(n,d,(re4(n/d)+I*im4(n/d))*d^2));
P4 = mk(n->sumdiv(n,d,(re4(n/d)-I*im4(n/d))*d^2));
print("R3 - (P3+P4)   = ", R3 - (P3+P4));
print("R4 - I*(P3-P4) = ", R4 - I*(P3-P4));
default(realprecision, 60);
s5=sqrt(5); ph5=(11+5*s5)/2; phm5=(5*s5-11)/2;
print("max|coeff| of  (R3 - phm5*R4) - ((1-I*phm5)*P3 + (1+I*phm5)*P4) = ", vecmax(vector(NQ-1,i,abs(polcoeff((R3-phm5*R4) - ((1-I*phm5)*P3+(1+I*phm5)*P4), i)))));
print("max|coeff| of  (R3 + ph5*R4)  - ((1+I*ph5)*P3 + (1-I*ph5)*P4)  = ", vecmax(vector(NQ-1,i,abs(polcoeff((R3+ph5*R4) - ((1+I*ph5)*P3+(1-I*ph5)*P4), i)))));
print("Phi_new  q-exp (as a+b*sqrt5): a = ", vector(9,i,polcoeff(R3,i-1)+11/2*polcoeff(R4,i-1)), "  b = ", vector(9,i,5/2*polcoeff(R4,i-1)));
quit;
