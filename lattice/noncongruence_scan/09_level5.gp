/* The level-5 Fricke non-congruence row (the Gamma_0(5)+5 analogue of
   Beukers 1987 Thm 3):
     u = (eta_5/eta_1)^6,  t = u/(1+22u+125u^2)  (Hauptmodul of Gamma_0(5)+5),
     F = E_{2,5} = (5E_2(5tau)-E_2(tau))/4 in M_2(Gamma_0(5)), F|_2 W_5 = -F,
     g = sqrt(F) weight one on an index-2 subgroup,  a_n = 2^n [t^n] g.        */
{
MQ=210; NROW=200;
E=buildE([1,5],MQ);
u=etaq(E,[-6,6],1,MQ);
F5 = 1 + 6*sum(n=1,MQ-1, (sigma(n) - if(n%5==0,5*sigma(n/5),0))*'q^n) + O('q^MQ);
tq = u/(1+22*u+125*u^2)+O('q^MQ);
print("t = ", Ser(vector(9,i,polcoeff(tq,i-1)),'q));
qt = serreverse(tq);
print("q(t) integral: ", vecmax(vector(NROW,i,denominator(polcoeff(qt,i))))==1);
Ft = subst(F5,'q,qt); g = sqrt(Ft);
c  = vector(NROW+1,i,polcoeff(g,i-1));
print("v_2(den c_n), n=1..12: ", vector(12,i,valuation(denominator(c[i+1]),2)));
a1 = vector(NROW+1,i,c[i]);            \\ lambda = 1
a2 = vector(NROW+1,i,2^(i-1)*c[i]);    \\ lambda = 2
print("lambda=1 integral? ", vecmax(vector(NROW+1,i,denominator(a1[i])))==1);
print("lambda=2 integral to n=",NROW,"? ", vecmax(vector(NROW+1,i,denominator(a2[i])))==1);
print("a_n = ", vector(8,i,a2[i]));
r = fitmin(a2,3,3);
print("minimal recurrence: order ",r[1]," degree ",r[2]);
print("  coeffs: ", r[3]~);
cp = charpol_from(r[3],r[1],r[2]);
print("  char poly ",cp,"  roots ",polroots(cp));
\\ companion and k
kk = compk(r,120);
print("k = ",kk[1]);
my(b=kk[2]);
print("  d_n^1 b_n integral? ", vecmax(vector(121,i,denominator(DL[i]^1*b[i])))==1);
print("  d_n^2 b_n integral? ", vecmax(vector(121,i,denominator(DL[i]^2*b[i])))==1);
\\ Casoratian
print("  W_n (n=0..4): ", vector(5,i,a2[i]*b[i+1]-a2[i+1]*b[i]));
}
