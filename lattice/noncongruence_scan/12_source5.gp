/* Theorem R3 check for the level-5 Fricke row:
   Psi = g^3 t / lambda  (weight three),  Theta = theta_q^{-2} Psi,
   b_n = lambda^n [t^n] (g * Theta)  should equal the census companion.        */
{
default(parisize,3000000000);
MQ=60; NR=40; LAM=2;
E=buildE([1,5],MQ);
u=etaq(E,[-6,6],1,MQ);
F5 = 1 + 6*sum(n=1,MQ-1, (sigma(n) - if(n%5==0,5*sigma(n/5),0))*'q^n) + O('q^MQ);
tq = u/(1+22*u+125*u^2)+O('q^MQ);
f  = sqrt(F5);                       \\ weight one, as a q-series
Psi = f^3*tq/LAM;
print("Psi = ", Ser(vector(7,i,polcoeff(Psi,i-1)),'q));
print("lam^m psi(m) integral for m<=", MQ-1, ": ",
      vecmax(vector(MQ-1,m,denominator(LAM^m*polcoeff(Psi,m))))==1);
Theta = sum(m=1,MQ-1, polcoeff(Psi,m)/m^2*'q^m) + O('q^MQ);
qt = serreverse(tq);
B  = subst(f*Theta,'q,qt);
bR = vector(NR+1); bR[1]=0; bR[2]=1;
for(n=1,NR-1, my(P=88*n^2+44*n+6, Q=-64*n^2+64*n-12);
  bR[n+2]=(P*bR[n+1]-Q*bR[n])/(n+1)^2);
bS = vector(NR+1,i, LAM^(i-1)*polcoeff(B,i-1));
print("b from the source  : ", vector(6,i,bS[i]));
print("b from the recurrence: ", vector(6,i,bR[i]));
print("agree for n<=",NR,": ", bS==bR);
}
