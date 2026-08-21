\\ rows_07_modular.gp -- the modular identity  sum a_n t^n = theta(q)  for the three
\\ lambda=1 rows (Domb, T, Cooper s7).  t is recovered by exact series reversion from
\\ the Sym^1 row and the weight-one CM theta series, then tested against the candidate
\\ eta quotient / identified by an exact linear solve on log(t/q).
default(parisizemax,6000000000);
M=201;
E(d) = prod(n=1,(M-1)\d, 1 - q^(d*n)) + O(q^M);

srow(c1,c2,c3,d1,d2,d3,B,N) = { my(a=vector(N+2)); a[1]=1; a[2]=B;
  for(n=1,N, a[n+2]=((c1*n^2+c2*n+c3)*a[n+1]-(d1*n^2+d2*n+d3)*a[n])/(n+1)^2); a; }

thD = E(1)^2*E(3)^2/(E(2)*E(6));                                    \\ level 12, chi_-3
thT = 1 + 2*sum(n=1,M-1, sumdiv(n,d,kronecker(-8,d))*q^n) + O(q^M); \\ theta_{x^2+2y^2}
th7 = 1 + 2*sum(n=1,M-1, sumdiv(n,d,kronecker(-7,d))*q^n) + O(q^M); \\ theta_{x^2+xy+2y^2}

tD = -q*(E(2)*E(6)/(E(1)*E(3)))^6 + O(q^M);       \\ Domb curve parameter, sign-flipped
tT =  q*(E(1)*E(8)/(E(2)*E(4)))^8 + O(q^M);       \\ T curve parameter (SPORADIC_SCAN2 sec.4)

revert(a, N, th) = { my(A = sum(n=0,N,a[n+1]*'x^n) + O('x^(N+1)));
  subst(serreverse(A-1), 'x, th - 1); }

\\ exact eta-quotient test on the divisor list DS:  log(t/q) = sum_d r_d log(E_d),
\\ sum_d r_d*d = 24.  Solve on the first #DS coefficient equations, then verify.
etaid(t, DS) = { my(nd=#DS, K=nd+8, LH=log(t/q+O(q^(K+2))), Ms, v, r);
  Ms = matrix(K+1, nd);  v = vector(K+1)~;
  for(j=1,nd, my(l=log(E(DS[j])+O(q^(K+2)))); for(i=1,K, Ms[i,j]=polcoef(l,i)));
  for(i=1,K, v[i]=polcoef(LH,i));
  for(j=1,nd, Ms[K+1,j]=DS[j]); v[K+1]=24;
  r = matinverseimage(Ms, v);
  if(type(r)!="t_COL" || #r==0, return([0, 0]));
  [r~, Ms*r == v]; }

report(nm, a, th, tcand, DSlist, N) = { my(t=revert(a,N,th), P=180);
  print("=== ", nm, " ===");
  print("  theta(q) = ", th + O(q^11));
  print("  t(q)     = ", t + O(q^11));
  print("  t in Z[[q]] up to q^",P,"? ", vector(P,i,denominator(polcoef(t,i)))==vector(P,i,1));
  if(tcand!=0,
     print("  t equals the stated eta quotient up to q^",P,"? ",
           vector(P,i,polcoef(t,i))==vector(P,i,polcoef(tcand,i))));
  for(j=1,#DSlist, my(r=etaid(t,DSlist[j]));
     print("  eta-quotient fit on divisors ",DSlist[j],": exponents ",r[1],
           "   overdetermined system consistent: ",r[2]));
  print("  log|c_n|/n for n=20,50,100,150,180: ",
        vector(5,i,my(n=[20,50,100,150,180][i]); log(abs(polcoef(t,n)))/n));
  print("");
}

aD = srow(20,10,2, 64,-64,16, 2, M+5);
aT = srow(24,12,2, 16,-16,4, 2, M+5);
a7 = srow(26,13,2, -27,27,-6, 2, M+5);

default(realprecision,40);
report("Domb  (level 12, disc -3)", aD, thD, tD, [[1,2,3,6],[1,2,3,4,6,12]], M-4);
report("T     (level 8,  disc -8)", aT, thT, tT, [[1,2,4,8]],               M-4);
report("s7    (disc -7)",           a7, th7, 0,  [[1,7],[1,2,7,14],[1,3,7,21],[1,2,3,6,7,14,21,42]], M-4);

\\ ---- s7: t is NOT an eta quotient, but IS a degree-2 rational function of the
\\ Gamma_0(7) hauptmodul h = (eta_1/eta_7)^4.  Fit t*Q(h) = P(h), deg <= D.
{
my(N=M-4, A=sum(n=0,N,a7[n+1]*'x^n)+O('x^(N+1)), t, h=(E(1)/E(7))^4/q, D, K, cols, Ms, ker, v, P, Q, tt);
t = subst(serreverse(A-1),'x, th7-1);
for(D=1,4,
  K=4*D+40; cols=List();
  for(i=0,D, listput(cols, vector(K,m, polcoef(h^i+O(q^K), m-1-D))));
  for(j=0,D, listput(cols, vector(K,m, polcoef(-t*h^j+O(q^K), m-1-D))));
  Ms=Mat(vector(#cols,c,Vec(cols[c])~)); ker=matker(Ms);
  print("  s7: t = P(h)/Q(h) with deg <= ",D,"?  dim ker = ",#ker));
D=2; K=4*D+40; cols=List();
for(i=0,D, listput(cols, vector(K,m, polcoef(h^i+O(q^K), m-1-D))));
for(j=0,D, listput(cols, vector(K,m, polcoef(-t*h^j+O(q^K), m-1-D))));
Ms=Mat(vector(#cols,c,Vec(cols[c])~)); ker=matker(Ms); v=ker[,1]; v=v/content(v);
P=sum(i=0,D,v[i+1]*'X^i); Q=sum(j=0,D,v[D+2+j]*'X^j);
print("  s7: P(X) = ",P,"   Q(X) = ",Q);
tt = subst(P,'X,h)/subst(Q,'X,h);
print("  s7: t == P(h)/Q(h) verified to q^120? ",
      vector(120,i,polcoef(t,i))==vector(120,i,polcoef(tt,i)));
print("  s7: equivalently with u = 1/h = (eta_7/eta_1)^4 in q+q^2 Z[[q]]:");
print("        t = u/(1 + 13u + 49u^2)   -- so t in q + q^2 Z[[q]], hence a_n in Z.");
}
\q
