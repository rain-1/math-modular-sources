/* dwork.gp -- Dwork/Gauss congruence test  a_{m p^s} = a_{m p^{s-1}} mod p^s
   (the h=1 shadow of Beukers-Vlasenko I Cor 5.9 / Mellit-Vlasenko Thm 1).
   Tested on: the census rows (integral), Zudilin's Catalan row rescaled
   R_m = 16^m Q_m, and the Brown-Zudilin zeta(5) row Q_n.  */
default(parisizemax, 6000000000);
read("lattice/euler_criterion/rows.gp");
N = 3000;

/* Brown-Zudilin row (arXiv:2210.03391 sect.2, symmetric case) */
c0(n) = 2*(2*n+1)*(41218*n^3-48459*n^2+20010*n-2871)*(n+1)^5;
c1(n) = -(97604224*n^9 + 178061760*n^8 + 72005308*n^7 - 48634688*n^6 - 39076836*n^5 + 2622730*n^4 + 7581006*n^3 + 920112*n^2 - 543402*n - 120582);
c2(n) = -2*n*(3874492*n^8 - 2617900*n^7 - 3144314*n^6 + 2947148*n^5 + 647130*n^4 - 1182926*n^3 + 115771*n^2 + 170716*n - 44541);
c3(n) = n*(41218*n^3+75195*n^2+46746*n+9898)*(n-1)^5;
{ bzrun(v0,v1,v2,M) = my(v=vector(M+1)); v[1]=v0; v[2]=v1; v[3]=v2;
  for(n=2,M-1, v[n+2] = -(c1(n)*v[n+1] + c2(n)*v[n] + c3(n)*v[n-1])/c0(n)); v; }

test(nm, A, M) = {
  my(bad=0, tot=0, firstbad=0, nonint=0);
  for(i=1,#A, if(denominator(A[i])!=1, nonint++));
  for(j=1,4, my(p=[2,3,5,7][j]);
    my(s=1); while(p^s<=M,
      for(m=1, M\p^s, my(x=A[m*p^s+1]-A[m*p^(s-1)+1]);
        tot++;
        if(x!=0 && valuation(x,p)<s, bad++; if(firstbad==0, firstbad=[p,s,m])));
      s++));
  print(nm, ": non-integral entries=", nonint, "   Dwork tests=", tot, "  failures=", bad,
        if(bad>0, Str("  first: (p,s,m)=", firstbad), ""));
}

R = rowR2(10,3,9,N); test("row C  a_n            ", R[1], N);
R = rowR2(12,4,32,N); test("row E  a_n            ", R[1], N);
R = rowR2(17,6,72,N); test("row F  a_n            ", R[1], N);
R = rowR3(17,5,1,N);  test("Apery zeta(3) a_n     ", R[1], N);
R = rowR3(10,4,64,N); test("Domb alpha a_n        ", R[1], N);
R = rowS18(N);        test("Cooper s18 a_n        ", R[1], N);
R = rowCusp(N);       test("cusp row a_n          ", R[1], N);

MZ = 1200;
Z = rowZud(MZ);
{ my(Rm=vector(MZ+1,i, 16^(i-1)*Z[1][i]));
  test("Zudilin 16^m Q_m      ", Rm, MZ);
  print("  v_2(16^m Q_m) for m=1..12: ", vector(12,i,valuation(Rm[i+1],2)));
  print("  s_2(m)*2 for m=1..12     : ", vector(12,i,2*hammingweight(i)));
  print("  R_0..R_6 = ", vector(7,i,Rm[i])); }

MB = 900;
{ my(Q=bzrun(1,21,2989,MB));
  test("Brown-Zudilin Q_n     ", Q, MB); }
