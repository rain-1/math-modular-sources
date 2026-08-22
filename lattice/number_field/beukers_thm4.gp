default(parisizemax, 4000000000);
default(realprecision, 120);
\\ ---------------------------------------------------------------
\\ Task 1: Beukers, Irrationality proofs using modular forms,
\\ Asterisque 147-148 (1987), Theorem 4 -- exact reproduction.
\\ K = Q(sqrt5), w = (1+sqrt5)/2 realised as Mod(y, y^2-y-1).
\\ ---------------------------------------------------------------
w   = Mod(y, y^2-y-1);
s5  = 2*w-1;                  \\ sqrt5
print("check s5^2 = ", lift(s5^2));

\\ coefficient extraction in the basis [1,w]
co(e) = my(p=lift(e)); [polcoeff(p,0), polcoeff(p,1)];
isint(e) = my(v=co(e)); (denominator(v[1])==1) && (denominator(v[2])==1);

\\ embeddings
ev1(e) = my(v=co(e)); v[1] + v[2]*(1+sqrt(5))/2;
ev2(e) = my(v=co(e)); v[1] + v[2]*(1-sqrt(5))/2;

\\ ---------------- (a) the recurrence and printed values -------------
A = 124+55*s5;
B = 34+15*s5;
Pn(n) = (A*n*(n+1) + B)*(2*n+1);
\\ (n+1)^3 d_{n+1} = Pn(n) d_n - n^3 d_{n-1}
NMAX = 320;
dd = vector(NMAX+2);
dd[1] = 1;                    \\ dd[i] = d_{i-1}
\\ n=0 : 1*d_1 = B*d_0
dd[2] = Pn(0)*dd[1];
for(n=1, NMAX, dd[n+2] = (Pn(n)*dd[n+1] - n^3*dd[n])/(n+1)^3);
print("d_0 = ", lift(dd[1]), "   (in basis 1,w)");
for(j=1,3, my(v=co(dd[j+1])); \\ convert a+b*w  ->  A + B*sqrt5 : a+b*w = (a+b/2) + (b/2)sqrt5
   print("d_",j," = ", v[1]+v[2]/2, " + ", v[2]/2, "*sqrt5"));
print("printed d_1 = 34 + 15 sqrt5 : ", dd[2]==34+15*s5);
print("printed d_2 = 7111 + 3180 sqrt5 : ", dd[3]==7111+3180*s5);
print("printed d_3 = 2040334 + 912465 sqrt5 : ", dd[4]==2040334+912465*s5);

\\ ---------------- (b) integrality of d_n ---------------------------
bad = -1;
for(n=0, 300, if(!isint(dd[n+1]), bad=n; break));
print("first non-integral d_n for n<=300 : ", if(bad<0,"NONE (all integral in Z[w])",bad));
\\ also: are they in Z[sqrt5] (i.e. b even)?
bad2 = -1;
for(n=0, 300, my(v=co(dd[n+1])); if(v[2]%2!=0, bad2=n; break));
print("first d_n NOT in Z[sqrt5] (odd w-coefficient), n<=300 : ", if(bad2<0,"NONE",bad2));

\\ ---------------- (c) characteristic polynomial ---------------------
\\ derive: divide recurrence by n^3 and let n->oo:  lam^2 - 2A lam + 1 = 0
print("derived char poly : lam^2 - (", lift(2*A), ") lam + 1   [2*A = ", 2*(124+55*2.23606797749978859), "]");
print("  2A in a+b sqrt5 form: ", 248, " + ", 110, "*sqrt5   ; matches 248+110sqrt5 : ", 2*A == 248+110*s5);
default(realprecision, 60);
tr1 = ev1(2*A); tr2 = ev2(2*A);
r1a = (tr1+sqrt(tr1^2-4))/2; r1b = (tr1-sqrt(tr1^2-4))/2;
r2a = (tr2+sqrt(tr2^2-4))/2; r2b = (tr2-sqrt(tr2^2-4))/2;
print("v1 (sqrt5->+2.2360679...) trace = ", tr1);
print("  roots: ", r1a, "  and  ", r1b);
print("v2 (sqrt5->-2.2360679...) trace = ", tr2);
print("  roots: ", r2a, "  and  ", r2b);
print("product of roots at v1 = ", r1a*r1b, " ; at v2 = ", r2a*r2b);

\\ ---------------- (d) Beukers' accounting ---------------------------
rho1 = r1a;          \\ fold at v1: 1/|lambda_small| = lambda_big
rho2 = r2b;          \\ no fold at v2: 1/|lambda_big| = lambda_small
gm  = sqrt(rho1*rho2);
print("rho_v1 = ", rho1);
print("rho_v2 = ", rho2, "   (= 1/", r2a, ")");
print("geometric mean sqrt(rho1*rho2) = ", gm);
print("e^3 = ", exp(3));
print("margin (log) = (1/2)(log rho1 + log rho2) - 3 = ", log(gm)-3);
print("ratio  gm/e^3 = ", gm/exp(3));

\\ ---------------- (e) counterfactual both-folded --------------------
rho2p = r2a;
gmp = sqrt(rho1*rho2p);
print("counterfactual rho_v2' = ", rho2p);
print("  gm' = sqrt(", rho1*rho2p, ") = ", gmp);
print("  margin' = ", log(gmp)-3, "   ratio' = ", gmp/exp(3));
\\ lambda_2 as an element of K : it is a root of lam^2-2A lam+1 in K?
disc = (2*A)^2-4;
print("disc of char poly = ", lift(disc), " ; is it a square in K? ", issquare(disc));
