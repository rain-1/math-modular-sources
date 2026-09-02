\\ 30_calib.gp -- Task (1): calibration of the Shimura-Borcherds dictionary
\\ for Pasol-Zudilin's level-one magnetic forms F4a, F4b (k=2, D=1)
\\ and the doubly magnetic F6 (k=3, D=-3).   All arithmetic exact.
default(parisize, 2000000000);
NN = 320;          \\ series precision
MM = 300;          \\ range of m

qz = 'q + O('q^NN);
{ E4s = 1 + 240*sum(n=1, NN-1, sigma(n,3)*'q^n) + O('q^NN); }
{ E6s = 1 - 504*sum(n=1, NN-1, sigma(n,5)*'q^n) + O('q^NN); }
DEL = (E4s^3 - E6s^2)/1728;

F4a = DEL/E4s^2;
F4b = E4s*DEL/E6s^2;
F6  = E6s*DEL/E4s^3;

print("== sanity: leading coefficients ==");
print("F4a q^1..q^5 : ", vector(5,m,polcoeff(F4a,m)));
print("F4b q^1..q^5 : ", vector(5,m,polcoeff(F4b,m)));
print("F6  q^1..q^5 : ", vector(5,m,polcoeff(F6 ,m)));

\\ ---------- k = 2 : A(n) = sum_{d|n} d a(n^2/d^2),  beta = (A/n) * mu ----------
{ Avec(F,M) = vector(M, m, polcoeff(F,m)); }
{ Betak2(A) = my(M=#A, ap, b);
  ap = vector(M, n, A[n]/n);
  b  = vector(M);
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*ap[n/d]); b[n]=s);
  [ap, b];
}
\\ ---------- k = 3 : A(n) = sum_{d|n} (-3/d) d^2 a(3n^2/d^2), beta3 = (A/n^2)*(mu.chi) ----
{ Betak3(A) = my(M=#A, ap, b);
  ap = vector(M, n, A[n]/n^2);
  b  = vector(M);
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*kronecker(-3,d)*ap[n/d]); b[n]=s);
  [ap, b];
}

{ report(nam, A, ap, b, k) =
  my(M=#A, badA=0, badb=0, sq, worst, tv, dist, g);
  print("");
  print("=========== ", nam, "  (k=", k, ") ===========");
  \\ magnetism: n^k' | A(n)  (k'=1 for k=2, k'=2 for k=3)
  for(n=1,M, if(denominator(ap[n])!=1, badA=n; break));
  print("  A(n)/n^", if(k==2,1,2), " integral for all n<=", M, ": ", if(badA==0,"YES",concat("FAIL at n=",badA)));
  for(n=1,M, if(denominator(b[n])!=1, badb=n; break));
  print("  beta(n) integral for all n<=", M, ": ", if(badb==0,"YES",concat("FAIL at n=",badb)));
  print("  beta(1..14) = ", vector(14,i,b[i]));
  \\ a(|D| m^2) = m^{k-1} beta(m)
  sq = vector(M, m, m^(k-1)*b[m]);
  print("  a(|D|*m^2) for m=1..8 : ", vector(8,i,sq[i]));
  \\ exact power of m dividing a(|D| m^2)
  tv = vector(M);
  for(m=2,M, my(t=0, x=sq[m]); if(x==0, t=-1, while(x%(m^(t+1))==0, t++)); tv[m]=t);
  dist = vector(6);
  for(m=2,M, my(t=tv[m]); if(t>=0 && t<=5, dist[t+1]++));
  print("  #\{m in [2,",M,"] : m^t || a(|D|m^2)\}  t=0,1,2,3,4,5 : ", dist);
  print("  #\{m : a(|D|m^2)=0\} = ", sum(m=2,M, tv[m]==-1));
  print("  m with t>=2 : ", select(x->x, vector(M,m, if(m>1 && tv[m]>=2, m, 0))));
  \\ gcd of beta(m) over m in a window
  print("  gcd_\{2<=m<=",M,"\} beta(m) = ", gcd(vector(M-1,i,b[i+1])));
  print("  gcd_\{2<=m<=",M,"\} a(|D|m^2)/m^",k-1," = ", gcd(vector(M-1,i,b[i+1])));
}

A1 = Avec(F4a,MM); T1 = Betak2(A1); report("F4a = Delta/E4^2", A1, T1[1], T1[2], 2);
A2 = Avec(F4b,MM); T2 = Betak2(A2); report("F4b = E4 Delta/E6^2", A2, T2[1], T2[2], 2);
A3 = Avec(F6 ,MM); T3 = Betak3(A3); report("F6 = E6 Delta/E4^3", A3, T3[1], T3[2], 3);

print("");
print("== square-index coefficients a(|D| m^2) = m^(k-1) beta(m), m=1..8 ==");
{ for(i=1,8, print("  m=",i,
    "   F4a a(m^2)=", i*T1[2][i],
    "   F4b a(m^2)=", i*T2[2][i],
    "   F6  a(3m^2)=", i^2*T3[2][i])); }

print("");
print("== v_p(a(p^{2k})) for small p,k ; a(p^{2k}) = p^k beta(p^k) [k=2 rows], p^{2k} beta(p^k) [F6] ==");
{ for(ip=1,6, my(p=prime(ip));
    for(e=1, floor(log(MM)/log(p)),
      my(m=p^e);
      print("  p=",p," e=",e,
        "  v_p(beta_F4a(p^e))=", if(T1[2][m]==0,"inf",valuation(T1[2][m],p)),
        "  v_p(beta_F4b(p^e))=", if(T2[2][m]==0,"inf",valuation(T2[2][m],p)),
        "  v_p(beta_F6(p^e))=",  if(T3[2][m]==0,"inf",valuation(T3[2][m],p)))));
}

write("30_beta_F4a.txt", T1[2]);
write("30_beta_F4b.txt", T2[2]);
write("30_beta_F6.txt",  T3[2]);
write("30_A_F4a.txt", A1);
write("30_A_F4b.txt", A2);
write("30_A_F6.txt",  A3);
quit;
