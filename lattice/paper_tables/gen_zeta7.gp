\\ Level-24 zeta(7) Apery system, from
\\ "chatgpt-research-archive/Closed binomial-harmonic forms for a level-24 zeta(7) Apery system.md"
\\ and eta-quotient r(tau) from modular_apery_cm_isogeny_research_notes.txt sec VII.C.
\\
\\ A(z) = sqrt(1-34s+s^2) * Apery(s)^3,  s=(z/(1-z))^2   [A_0=1, integral]
\\ r(tau) = [eta(2t)eta(12t)/(eta(4t)eta(6t))]^6 = q*(prod ratio)^6,  z=r/(1+r)
\\ Phi(q)=sum g_m q^m,  g_m=sum_{d|24,d|m} c_d sigma_7(m/d),
\\   (c_1,c_2,c_3,c_4,c_6,c_8,c_12,c_24)=(1,-588,11583,-27456,-138996,585728,-762048,331776)
\\ Psi(q)=sum g_m/m^7 q^m = D^{-7}Phi
\\ B(z(q)) = A(z(q)) Psi(q)   =>   B(z) = A(z) * Psi(q(z)),  q(z)=serreverse(z(q))
default(realprecision,60);
N=220; \\ series order in q and z

\\ ---- 1. r(q), z(q), q(z) from eta quotients ----
fk(k,NN)={my(r=1+O(q^(NN+1))); for(n=1, NN\k+1, r*=(1-q^(k*n))); r};
rq = q*(fk(2,N)*fk(12,N)/(fk(4,N)*fk(6,N)))^6;
zq = rq/(1+rq);
Qz = serreverse(zq); \\ q as a power series in z (variable still called q)

\\ ---- 2. Phi(q), Psi(q) ----
cd = [1,-588,11583,-27456,-138996,585728,-762048,331776];
dl = [1,2,3,4,6,8,12,24];
gvec = vector(N, m, sum(i=1,#dl, if(m%dl[i]==0, cd[i]*sigma(m\dl[i],7), 0)));
Phi = sum(m=1,N, gvec[m]*q^m) + O(q^(N+1));
Psi = sum(m=1,N, (gvec[m]/m^7)*q^m) + O(q^(N+1));

\\ ---- 3. A(z) via Apery(s)^3 * sqrt(1-34s+s^2), s=(z/(1-z))^2 ----
Ms = N\2+3;
apery = vector(Ms+1, mp1, my(m=mp1-1); sum(k=0,m, binomial(m,k)^2*binomial(m+k,k)^2));
Aseries_s = sum(m=0,Ms, apery[m+1]*q^m) + O(q^(Ms+1)); \\ reuse var q as s here (careful: separate poly)
\\ use a fresh variable for s to avoid clashing with q used elsewhere
Aseries_s = subst(Aseries_s, q, 's);
sqrtpart = sqrt(1 - 34*'s + 's^2 + O('s^(Ms+1)));
Acubed_s = Aseries_s^3 * sqrtpart; \\ series in s to O(s^(Ms+1))

\\ substitute s=(z/(1-z))^2 as a series in z (z itself represented via variable q, since Qz etc use q)
zser = q + O(q^(N+1));
sser = (zser/(1-zser))^2;
Aser = subst(Acubed_s, 's, sser); \\ A(z) as series in q(=z), truncated automatically to O(q^(N+1))

\\ ---- 4. B(z) = A(z) * Psi(q(z)) ----
PsiOfQz = subst(Psi, q, Qz);
Bser = Aser * PsiOfQz;

NN = N-2; \\ safe margin below series precision
An = vector(NN+1, i, polcoeff(Aser, i-1));
Bn = vector(NN+1, i, polcoeff(Bser, i-1));

print("A_0..A_10 = ", vector(11,i,An[i]));
print("B_1..B_10 = ", vector(10,i,Bn[i+1]));
print("series order available: NN=",NN);

\\ ---- 5. Cross-check A_n against the boxed binomial-harmonic formula ----
\\ y_l = -1_{l=1} - 2 sum_{p=1}^l 8^p/p C(2p-2,p-1) C(l+p-2,l-p)
yl(l)={if(l==0, return(1)); \\ y_0 from sqrt(1)=1 series matches apery convention below separately
 -( (l==1) + 2*sum(p=1,l, 8^p/p*binomial(2*p-2,p-1)*binomial(l+p-2,l-p)))};
Msmall=8;
ylist = vector(Msmall+1, i, yl(i-1));
alist = vector(Msmall+1, i, apery[i]);
Cm(m)={sum(l=0,m, sum(i=0,m-l, sum(j=0,m-l-i, my(k=m-l-i-j); if(k>=0, ylist[l+1]*alist[i+1]*alist[j+1]*alist[k+1], 0))))};
Anboxed(n)={if(n==0, return(1)); sum(m=1,n\2, binomial(n-1,2*m-1)*Cm(m))};
print("boxed A_n check, n=0..8: ", vector(9,i,Anboxed(i-1)));
print("z-route     A_n,  n=0..8: ", vector(9,i,An[i]));

\\ ---- 6. integrality: d_n^7 B_n in Z ? ----
{my(ok=1, dn=1); for(n=1,NN, dn=lcm(dn,n); my(val=dn^7*Bn[n+1]);
  if(val!=round(val), ok=0; print("FAIL n=",n," val=",val)));
 print("d_n^7*B_n integral for all n=1..",NN,"? ",ok);}

\\ ---- 7. limit B_n/A_n -> (1463/13824) zeta(7) ----
target = (1463/13824)*zeta(7);
lim = 1.0*Bn[NN+1]/An[NN+1];
print("B_n/A_n at n=",NN," = ",lim);
print("target (1463/13824)zeta(7) = ",1.0*target);
print("diff = ",lim-target);
print("|diff|^(1/n) ~ ",abs(lim-target)^(1.0/NN)," (expect 1/sqrt(2)=",1/sqrt(2.),")");

\\ ---- 8. p-adic slopes of B_n/A_n and v_p(A_n) growth, p=2,3,5,7 ----
{my(ns=[30,60,90,120,150,180,NN]);
 for(pi=1,4, my(p=[2,3,5,7][pi]);
   my(diffs=vector(#ns-1,k, my(n=ns[k+1]); valuation(Bn[n+1]/An[n+1]-Bn[n]/An[n],p)));
   print("p=",p," v_p(B_n/A_n - B_{n-1}/A_{n-1}) at n=",vecextract(ns,"2.."),": ",diffs);
   my(mv=vector(#ns,k,my(n=ns[k]); valuation(An[n+1],p)));
   print("     v_p(A_n) at n=",ns,": ",mv));}

\\ ---- 9. characteristic-root discussion (analytic, per coordinator's note) ----
lam1 = 4+2*sqrt(2.);
lam2 = lam1/sqrt(2.);
print("\nassuming dominant root lambda1=4+2sqrt2=",lam1,", and archimedean decay rate 1/sqrt2 = lambda2/lambda1,");
print("second-dominant root lambda2 = lambda1/sqrt2 = ",lam2," = 2+2sqrt2 numerically check: ",2+2*sqrt(2.));
\q
