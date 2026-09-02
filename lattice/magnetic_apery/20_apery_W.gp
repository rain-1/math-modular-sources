default(parisizemax,12000000000);
read("lib.gp");
NQ = 405; MA=60; MB=150; MC=400;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
print("### Apery host N=6, C=72, B=17.  Systematic scan in W = 72u + 1/u.");
print("### criterion: the shortest vector of the lattice L(M)={alpha: Xi_alpha in Z[[q]] mod q^M}");
print("### must have the SAME norm at M=", MA, " and M=", MB, ", and then Xi must be integral to m=", MC, ".");
{
tryh(hW, mult, wide, lab) = my(Qp, qv, Xs, T, res);
  Qp = subst(touW(hW, C), U, x)^mult;
  qv = vector(poldegree(Qp)+1, j, polcoeff(Qp, j-1));
  if(wide, Xs = xibasisW(us, F2, qv, MC), Xs = xibasisU(us, F2, qv, MC));
  if(type(Xs)=="t_INT", return(0));
  T = xitomat(Xs, MC);
  res = scanT(T, MA, MB, MC);
  if(res[1]==1, print("  *** MAGNETIC *** ", lab, "  dim=", #Xs, "  alpha=", res[2]); return(1));
  0;
}
print("");
print("-- A1: single rational pole pair, Q = (W+b)^2, anti-invariant (dim 1), |b| <= 800");
{for(b=-800, 800, tryh(W+b, 2, 0, concat("h=W+",b)));}
print("   done A1");
print("-- A2: same, widest space (dim 4), |b| <= 800");
{for(b=-800, 800, tryh(W+b, 2, 1, concat("h=W+",b)));}
print("   done A2");
print("-- A3: pole orders 1 and 3, widest space, |b| <= 300");
{for(b=-300, 300, tryh(W+b, 1, 1, concat("h=W+",concat(b,", order 1"))); tryh(W+b, 3, 1, concat("h=W+",concat(b,", order 3"))));}
print("   done A3");
print("-- A4: non-monic single pole, Q = (aW+b)^2, widest space, 2<=a<=40, |b|<=40a, gcd=1");
{for(a=2, 40, for(b=-40*a, 40*a, if(gcd(a,b)!=1, next); tryh(a*W+b, 2, 1, concat(concat("h=",a),concat("W+",b)))));}
print("   done A4");
quit;
