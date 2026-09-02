/* 01_series.gp -- TASK 1.  c(m) for m<=400 from the eta quotients; two independent
   formulas for Phi; magnetism m|c(m); companion formula reproduces b_n.          */
default(parisize, 1000000000);
read("lib.gp");
M = 400; N = M+6;
{ for(k=1,3, my(R=ROWS[k], S, SS, Ph, Ph2, cv, cpv, qx, Fx, av, RR, E, P, bb);
  S = Setup(k,N); Ph = S[4]; Ph2 = Phi2(k,N);
  print("=== ", R[1], " (N=", R[2], ", B=", R[3], ", C=", R[4], ") ===");
  print("  Phi = F*Dx  minus  x sqrt(P) F^2 = u(1-Cu^2)F^2/(1+Bu+Cu^2)^2 : ", Ph-Ph2);
  cv = vector(M, m, polcoeff(Ph,m));
  cpv = vector(M, m, cv[m]/m);
  print("  c(m) integral for m<=", M, " : ", #select(y->type(y)!="t_INT",cv)==0);
  print("  m | c(m) for m<=", M, "  : ", #select(m->cv[m]%m!=0, vector(M,m,m))==0, "  [MAGNETIC]");
  print("  c(m),  1<=m<=16 : ", vector(16,m,cv[m]));
  print("  c'(m), 1<=m<=16 : ", vector(16,m,cpv[m]));
  write(concat(R[1],"_c.txt"),  concat(concat("c(m), m=1..400, row ", R[1]), ":"));
  write(concat(R[1],"_c.txt"),  cv);
  write(concat(R[1],"_cp.txt"), concat(concat("c'(m)=c(m)/m, m=1..400, row ", R[1]), ":"));
  write(concat(R[1],"_cp.txt"), cpv);
  /* companion cross-check */
  SS = Setup(k,50); qx = serreverse(SS[3]); Fx = subst(SS[2],'q,qx);
  RR = genrow(k,45);
  av = vector(46, j, polcoeff(Fx,j-1));
  print("  a_n = [x^n] F  vs recurrence, n<=45 : ", av==vector(46,j,RR[1][j]));
  E = matrix(45,45); P = Fx;
  for(m=1,45, P = P*qx; for(n=m,45, E[n,m] = polcoeff(P,n)));
  print("  e_{n,m} = [x^n](F q^m) integral    : ", #select(y->type(y)!="t_INT"&&y!=0, concat(Vec(E)))==0);
  bb = vector(45, n, sum(m=1,n, cv[m]/m^3*E[n,m]));
  print("  companion b_n = sum c(m) m^-3 e_{n,m}  vs recurrence, n<=45 : ", bb==vector(45,n,RR[2][n+1]));
  print()); }
quit;
