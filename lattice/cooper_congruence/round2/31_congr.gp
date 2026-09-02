\\ 31_congr.gp -- Task (1)/(3): the twisted congruence (S) for the level-one PZ forms.
\\ Prediction derived from  f|T_{p^2} = psi(p) p f + p^{2k-1} c g_{m0 p^2}:
\\   k=2:  A(pm) = psi(p) p A(m) mod p^3,  psi = chi_{D0}, D0 = disc of the pole
\\   k=3:  A(pm) = p^2 A(m) mod p^5
default(parisize, 4000000000);
NN = 820; MM = 800;
{ E4s = 1 + 240*sum(n=1, NN-1, sigma(n,3)*'q^n) + O('q^NN); }
{ E6s = 1 - 504*sum(n=1, NN-1, sigma(n,5)*'q^n) + O('q^NN); }
DEL = (E4s^3 - E6s^2)/1728;
F4a = DEL/E4s^2;  F4b = E4s*DEL/E6s^2;  F6 = E6s*DEL/E4s^3;
A4a = vector(MM,m,polcoeff(F4a,m));
A4b = vector(MM,m,polcoeff(F4b,m));
A6  = vector(MM,m,polcoeff(F6 ,m));

\\ ---- (S) test:  c'(pm) = psi(p) c'(m) mod p^2   with c'(m)=A(m)/m  (k=2)
{ testS2(nam, A, D0) = my(M=#A, cp, msg);
  cp = vector(M, n, A[n]/n);
  print("");
  print("--- (S) for ", nam, " with psi = chi_", D0, " : c1(pm) = psi(p) c1(m) mod p^2 ---");
  forprime(p=2, 60,
    my(bad=0, cnt=0);
    for(m=1, M\p, cnt++; if((cp[p*m] - kronecker(D0,p)*cp[m]) % p^2 != 0, bad++));
    msg = if(bad==0, "PASS", Str("FAIL ", bad, "/", cnt));
    print("   p=", p, "  psi(p)=", kronecker(D0,p), "  tested m<=", M\p, " : ", msg));
}
\\ ---- k=3 :  c''(pm) = c''(m) mod p^3  with c''(m)=A(m)/m^2
{ testS3(nam, A) = my(M=#A, cp);
  cp = vector(M, n, A[n]/n^2);
  print("");
  print("--- (S3) for ", nam, " : c''(pm) = c''(m) mod p^3 ---");
  forprime(p=2, 60,
    my(bad=0, cnt=0);
    for(m=1, M\p, cnt++; if((cp[p*m] - cp[m]) % p^3 != 0, bad++));
    print("   p=", p, "  tested m<=", M\p, " : ", if(bad==0,"PASS",Str("FAIL ",bad,"/",cnt))));
}
testS2("F4a = Delta/E4^2", A4a, -3);
testS2("F4a (control, psi=1)", A4a, 1);
testS2("F4b = E4 Delta/E6^2", A4b, -4);
testS2("F4b (control, psi=1)", A4b, 1);
testS3("F6 = E6 Delta/E4^3", A6);
testS2("F6 (control k=2 reading)", A6, -3);

\\ ---- twisted Mobius inverse and rad(n)^2 | beta_psi(n)
{ radn(n) = my(f=factor(n)); prod(i=1,#f~, f[i,1]); }
{ betapsi(A, D0, e) = my(M=#A, cp, b);
  cp = vector(M, n, A[n]/n^e);
  b = vector(M);
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*kronecker(D0,d)*cp[n/d]); b[n]=s);
  b; }
{ radtest(nam, b, pw) = my(M=#b, bad=0, first=0);
  for(n=2,M, if(b[n]%(radn(n)^pw)!=0, bad++; if(first==0,first=n)));
  print("   ", nam, " : rad(n)^", pw, " | beta_psi(n) for 2<=n<=", M, " -> ", if(bad==0,"ALL PASS", Str("FAIL ",bad," values, first n=",first)));
}
print("");
print("--- rad(n)^2 | beta_psi(n)  (Cooper's condition (iii)) ---");
b4a = betapsi(A4a,-3,1); radtest("F4a, psi=chi_{-3}", b4a, 2);
b4a1= betapsi(A4a, 1,1); radtest("F4a, psi=1       ", b4a1,2);
b4b = betapsi(A4b,-4,1); radtest("F4b, psi=chi_{-4}", b4b, 2);
b4b1= betapsi(A4b, 1,1); radtest("F4b, psi=1       ", b4b1,2);
print("");
print("--- rad(n)^3 | beta_psi(n)  (weight-6 analogue, c''=A/n^2) ---");
b6  = betapsi(A6,  1,2); radtest("F6, psi=1        ", b6, 3);
b6c = betapsi(A6, -3,2); radtest("F6, psi=chi_{-3} ", b6c,3);
print("");
print("beta_psi(1..12) F4a (psi=chi_-3): ", vector(12,i,b4a[i]));
print("beta_psi(1..12) F4b (psi=chi_-4): ", vector(12,i,b4b[i]));
print("beta_psi(1..12) F6  (psi=1, /n^2): ", vector(12,i,b6[i]));
print("");
print("--- sharpness: is n^2 | beta_psi(n)? (Cooper has this) ---");
{ for(j=1,3, my(b=[b4a,b4b,b6][j], nm=["F4a","F4b","F6"][j], bad=0, first=0);
    for(n=2,MM, if(b[n]%(n^2)!=0, bad++; if(first==0,first=n)));
    print("   ", nm, " n^2|beta: ", if(bad==0,"ALL PASS",Str("FAIL ",bad," first n=",first))));
}
write("31_betapsi_F4a.txt", b4a);
write("31_betapsi_F4b.txt", b4b);
write("31_betapsi_F6.txt",  b6);
quit;
