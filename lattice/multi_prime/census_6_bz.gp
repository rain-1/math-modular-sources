default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/census_util.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/multi_prime/out/census_6_bz.log";
NN = 300;
d0(n) = 2*(2*n+1)*(41218*n^3-48459*n^2+20010*n-2871)*(n+1)^5;
d1(n) = -(97604224*n^9 + 178061760*n^8 + 72005308*n^7 - 48634688*n^6 - 39076836*n^5 + 2622730*n^4 + 7581006*n^3 + 920112*n^2 - 543402*n - 120582);
d2(n) = -2*n*(3874492*n^8 - 2617900*n^7 - 3144314*n^6 + 2947148*n^5 + 647130*n^4 - 1182926*n^3 + 115771*n^2 + 170716*n - 44541);
d3(n) = n*(41218*n^3+75195*n^2+46746*n+9898)*(n-1)^5;
runbz(v0,v1,v2) = { my(v = vector(NN+1)); v[1]=v0; v[2]=v1; v[3]=v2;
  for(n = 2, NN-1, v[n+2] = -(d1(n)*v[n+1] + d2(n)*v[n] + d3(n)*v[n-1])/d0(n)); v; };
qv = runbz(1, 21, 2989);
pv = runbz(0, 87/4, 1190161/384);
phv = runbz(0, 101/4, 344923/96);

logit(fn, Str("### Brown-Zudilin zeta(5) cellular row, N=", NN, " ###"));
{ my(ok = 1, k1); for(n = 0, 7,
    k1 = sum(a=0,n, sum(b=0,n, binomial(n+a,n)*binomial(n,a)^2*binomial(n+b,n)*binomial(n,b)^2*binomial(n+a+b,n)));
    if(qv[n+1] != k1, ok = 0));
  logit(fn, Str("  closed-form check Q_n = double binomial sum, n<=7: ", ok)); }
logit(fn, Str("  Q_0..Q_4 = ", vector(5, i, qv[i])));

ns = [50, 100, 200, 299];
report(fn, "BZ zeta(5) row  P/Q (zeta(5) direction)", qv, pv, primes(15), ns);
report(fn, "BZ zeta(5) row  Phat/Q (zeta(3) direction)", qv, phv, [2,3,5,7], ns);
logit(fn, "DONE");
quit;
