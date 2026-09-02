\\ 56_borcherds.gp (v2) -- TASK 4: is  P = prod_n (1-q^n)^{tau(n)}  a meromorphic modular
\\ form on Gamma_0(N)?  Test on the logarithmic derivative
\\   D log P = -sum_N (sum_{d|N} d tau(d)) q^N   (exact; needs tau(d) only for d <= N).
\\ If  P = c q^h (eta-quotient) R(u)  with R rational and poles only at the cusps and at
\\ the CM points g(u)=0, then D log P lies in
\\   span_Q { 1 } + span_Q { E_2(d tau) : d|N } + span_Q { F u^i / g(u)^j }
\\ since D log q^h = h, D log eta(d tau) = (d/24) E_2(d tau), D log R(u) = F u R'(u)/R(u).
\\ (The constant is essential: it carries the q^h prefactor.)
read("50_lib.gp");
TAUF = ["52_tau_s7.txt","53_tau_s10.txt","54_tau_s18.txt"];
{ E2s(d,M) = 1 - 24*sum(n=1, M\d, sigma(n)*'q^(d*n)) + O('q^(M+1)); }
{ dlogprod(e, M) = my(sg = vector(M, n, my(s=0); fordiv(n,d, s += d*e[d]); s));
  -sum(n=1,M, sg[n]*'q^n) + O('q^(M+1)); }
{ testspan(k, e, M, tag) =
  my(N=LEV[k], R=ROWS[k], Bc=R[3], Cc=R[4], S, u, F, g, L, BAS, NAMES, MAT, ker, dv, v);
  L = dlogprod(e, M);
  S = Setup(k, M+10); u = S[1]; F = S[2];
  g = 1 + Bc*u + Cc*u^2;
  BAS = List(); NAMES = List();
  listput(BAS, 1 + O('q^(M+1))); listput(NAMES, "1");
  fordiv(N, d, listput(BAS, E2s(d,M)); listput(NAMES, concat("E2(",concat(Str(d),"t)"))));
  for(j=0,2, for(i=-2,4,
     if(i==0 && j==0, next);
     listput(BAS, F*u^i/g^j);
     listput(NAMES, concat(concat("F u^",Str(i)),concat("/g^",Str(j))))));
  BAS = Vec(BAS); NAMES = Vec(NAMES);
  MAT = matrix(M+1, #BAS+1, r, c, if(c<=#BAS, polcoeff(BAS[c], r-1), -polcoeff(L, r-1)));
  ker = matker(MAT);
  print("  [",tag,"]  basis ",#BAS,", equations ",M+1,", dim ker = ", #ker, ", rank of basis alone = ", matrank(matrix(M+1,#BAS,r,c,polcoeff(BAS[c],r-1))));
  dv = 0;
  for(t=1,#ker, if(ker[,t][#BAS+1] != 0, dv = t; break));
  if(dv==0,
    print("     ==> NO relation:  D log P is NOT in the span."),
    v = ker[,dv]/ker[,dv][#BAS+1];
    print("     ==> RELATION:  D log P =");
    for(c=1,#BAS, if(v[c]!=0, print("            + (",v[c],") * ",NAMES[c]))));
}
print("== POSITIVE CONTROL: e(n) = 24 for all n, so P = eta^24/q, D log P = E2 - 1 ==");
{ for(k=1,3, testspan(k, vector(40,n,24), 40, concat("control ", NAM[k]))); }
print();
print("== POSITIVE CONTROL 2: e(n) = 4 for n=1..M with the s7 u:  prod(1-q^n)^4 = u/q * (eta7/eta)^0 ... ==");
print("   (u = q prod((1-q^7n)/(1-q^n))^4, so e(n) = -4 for 7 nmid n, +... ) skipped; control 1 suffices.");
print();
print("== THE TEST ==");
{ for(k=1,3, my(tau=read(TAUF[k])); testspan(k, tau, #tau, NAM[k])); }
quit;
