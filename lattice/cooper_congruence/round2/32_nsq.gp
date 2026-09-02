\\ 32_nsq.gp -- is the STRONG Cooper divisibility n^2 | beta_psi(n) also true for
\\ Pasol-Zudilin's level-one forms, away from the bad prime?
default(parisize, 4000000000);
NN = 820; MM = 800;
{ E4s = 1 + 240*sum(n=1, NN-1, sigma(n,3)*'q^n) + O('q^NN); }
{ E6s = 1 - 504*sum(n=1, NN-1, sigma(n,5)*'q^n) + O('q^NN); }
DEL = (E4s^3 - E6s^2)/1728;
F4a = DEL/E4s^2;  F4b = E4s*DEL/E6s^2;  F6 = E6s*DEL/E4s^3;
A4a = vector(MM,m,polcoeff(F4a,m));
A4b = vector(MM,m,polcoeff(F4b,m));
A6  = vector(MM,m,polcoeff(F6 ,m));
{ betapsi(A, D0, e) = my(M=#A, cp, b);
  cp = vector(M, n, A[n]/n^e);
  b = vector(M);
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*kronecker(D0,d)*cp[n/d]); b[n]=s);
  b; }
b4a = betapsi(A4a,-3,1);
b4b = betapsi(A4b,-4,1);
b6  = betapsi(A6,  1,2);
\\ per-prime valuation test:  v_p(beta(n)) >= w*v_p(n)  for every n
{ vtest(nam, b, w, plim) = my(M=#b);
  print("");
  print("--- ", nam, " :  v_p(beta(n)) >= ", w, "*v_p(n) ? ---");
  forprime(p=2, plim,
    my(bad=0, cnt=0, first=0, mn=99);
    for(n=p, M, my(e=valuation(n,p)); if(e>0, cnt++;
      my(vv=if(b[n]==0, 10^6, valuation(b[n],p)));
      if(vv < w*e, bad++; if(first==0,first=n));
      if(vv - w*e < mn, mn = vv-w*e)));
    print("   p=", p, " : ", if(bad==0, Str("PASS  (min excess ", mn, ")"), Str("FAIL ",bad,"/",cnt," first n=",first))));
}
vtest("F4a, psi=chi_{-3}, w=2", b4a, 2, 40);
vtest("F4b, psi=chi_{-4}, w=2", b4b, 2, 40);
vtest("F6,  psi=1, c''=A/n^2, w=3", b6, 3, 40);
print("");
print("--- same, w=1 (the PZ theorem level) ---");
vtest("F4a w=1", b4a, 1, 20);
vtest("F4b w=1", b4b, 1, 20);
vtest("F6  w=2", b6, 2, 20);
print("");
print("--- n^2 | beta restricted to n coprime to the bad prime ---");
{ chk(nam, b, w, badp) = my(M=#b, bad=0, first=0);
  for(n=2,M, if(n%badp!=0, my(r=1); fordiv(n,d,); if(b[n]%(n^w)!=0, bad++; if(first==0,first=n))));
  print("   ", nam, " : n^",w," | beta(n) for (n,",badp,")=1, n<=",M," -> ", if(bad==0,"ALL PASS",Str("FAIL ",bad," first n=",first)));
}
chk("F4a psi=chi_-3", b4a, 2, 2);
chk("F4b psi=chi_-4", b4b, 2, 3);
chk("F6  psi=1     ", b6,  3, 6);
quit;
