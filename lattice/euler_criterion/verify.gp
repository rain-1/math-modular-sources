/* ================================================================
   verify.gp -- census verification of the Euler-factor criterion
        xi_p = - Q(p-eval) * kappa_p ,
        kappa_p = 1/2 L_p(w+1, psi*om^{-w})  if the co-divisor character
                  phi of the depleted Eisenstein series is trivial,
        kappa_p = 0                          otherwise (phi != 1),
        kappa_p = 0                          for a cuspidal source.
   Run:  gp -q lattice/euler_criterion/verify.gp
   Set TGT (target p-adic digits) and PRC before reading, or use defaults.
   ================================================================ */
default(parisizemax, 8000000000);
read("lattice/euler_criterion/lp.gp");
read("lattice/euler_criterion/rows.gp");

if(type(TGT)!="t_INT", TGT = 900);
PRC = TGT + 60;
print("### target digits TGT=",TGT,"   L_p precision PRC=",PRC);
gettime();

k3_2    = Lp(3,triv ,2,PRC)/2;
k3_3    = Lp(3,triv ,3,PRC)/2;
k2_2    = Lp(2,triv ,2,PRC)/2;
k2_2_12 = Lp(2,chi12,2,PRC)/2;
k2_3    = Lp(2,triv ,3,PRC)/2;
k5_3    = Lp(5,triv ,3,PRC)/2;
print("L_p time ",gettime()," ms");
{
print("kappa_3^{w=1} = 1/2 L_3(2,1)     v_3=",valuation(k3_2,3));
print("kappa_3^{w=2} = 1/2 L_3(3,1)     v_3=",valuation(k3_3,3));
print("kappa_2^{w=1} = 1/2 L_2(2,1)     v_2=",valuation(k2_2,2));
print("kappa_2^{w=1,chi12}              v_2=",valuation(k2_2_12,2));
print("kappa_2^{w=2} = 1/2 L_2(3,1)     v_2=",valuation(k2_3,2));
print("kappa_5^{w=2} = 1/2 L_5(3,1)     v_5=",valuation(k5_3,5));
print("3*kappa_3^{w=1} mod 3^40 = ",lift(3*k3_2+O(3^40)),"   (CONJ_D: 8386265965554334030)");
print();
}

/* Kummer-sequence cross-checks of every kappa used, at exact rationals
   kappa_m = -(1-psi(p)p^{k-1}) B_{k,psi}/(2k), k -> -w along k = k0+(p-1)p^m */
kum(nm, p, chi, w, ms) =
{ my(k, km, ka);
  ka = if(w==1, if(p==3, k3_2, if(chi==chim4, k2_2, k2_2_12)), if(p==2,k2_3,if(p==3,k3_3,k5_3)));
  print1("  Kummer ",nm,":");
  for(i=1,#ms,
    my(m=ms[i], e=eulerphi(if(p==2,4,p))*p^m);
    k = -w + e*ceil((w+2)/e);          /* k = -w mod (p-1)p^m, k>=1 */
    km = -(1 - if(chi[1]%p==0, 0, chival(chi,p)*p^(k-1)))*Bchi(chi,k)/(2*k);
    print1("  [m=",m," k=",k," v=",valuation(km-ka,p),"]"));
  print();
}

/* --------- the census --------- */
rep(nm, p, sl, Qev, kp, gen) =
{ my(N, rw, A, B, x, xm1, cauchy, pred, vd);
  N = ceil((TGT+40)/sl);
  rw = gen(N); A=rw[1]; B=rw[2];
  x   = B[N+1]/A[N+1];  xm1 = B[N]/A[N];
  cauchy = valuation(x-xm1, p);
  pred = -Qev*kp;
  vd = if(pred===0, valuation(x,p), valuation(x-pred,p));
  printf("%-22s p=%d  N=%-5d  Q=%-8s  cauchy=%-6d  v_p(xi_N-pred)=%-6d  %s\n",
     nm, p, N, Str(Qev), cauchy, vd,
     if(vd >= min(cauchy,TGT) - 3*logint(N,p) - 5, "OK", "*** MISMATCH ***"));
}

nocv(nm,p,rw) =
{ my(A=rw[1],B=rw[2],N=#A-1);
  printf("%-22s p=%d  v_p(increment), n=%d..%d step 6: %s\n", nm, p, N-30, N,
     Str(vector(6, i, my(n=N-6*(6-i)); valuation(B[n+1]/A[n+1]-B[n]/A[n], p))));
}

tow(nm,p,rw,aa) =
{ my(A=rw[1],B=rw[2],N=#A-1);
  for(i=1,#aa, my(a=aa[i], v=List(), s=1);
    while(a*p^s<=N, listput(v, valuation(B[a*p^s+1]/A[a*p^s+1]-B[a*p^(s-1)+1]/A[a*p^(s-1)+1],p)); s++);
    printf("%-16s p=%d a=%d  v_p(xi_{a p^s}-xi_{a p^{s-1}}): %s\n", nm,p,a,Str(Vec(v))));
}

{
print("=== rows with Eisenstein source and Euler-factor divisibility ===");
rep("A   (7,2,-8)",   2, 3, 1,    0,       (N)->rowR2(7,2,-8,N));
rep("E   (12,4,32)",  2, 5, -1,   k2_2,    (N)->rowR2(12,4,32,N));
rep("F   (17,6,72)",  2, 3, -1,   k2_2_12, (N)->rowR2(17,6,72,N));
rep("F   (17,6,72)",  3, 2, -5/4, k3_2,    (N)->rowR2(17,6,72,N));
rep("B   (9,3,27)",   3, 3, -1,   k3_2,    (N)->rowR2(9,3,27,N));
rep("C   (10,3,9)",   3, 2, -1,   k3_2,    (N)->rowR2(10,3,9,N));
rep("alpha Domb",     2, 6, -2/3, k2_3,    (N)->rowR3(10,4,64,N));
rep("eps (12,4,16)",  2, 4, -1/2, k2_3,    (N)->rowR3(12,4,16,N));
rep("delta (7,3,81)", 3, 4, -1/2, k3_3,    (N)->rowR3(7,3,81,N));
rep("zeta (9,3,-27)", 3, 3, 1,    0,       (N)->rowR3(9,3,-27,N));
rep("eta (11,5,125)", 5, 3, -1,   k5_3,    (N)->rowR3(11,5,125,N));
print();
print("=== rows outside the twelve-family table ===");
rep("s18 Cooper",     3, 1, -1,   k3_2,    (N)->rowS18(N));
rep("Zudilin Catalan",2, 8, -2,   k2_2,    (N)->rowZud(N));
print("=== cuspidal sources: predicted xi_p = 0 ===");
rep("cusp L(f,2) L12",2, 2, 1,    0,       (N)->rowCusp(N));
print();
print("=== Kummer-sequence cross-checks of kappa (exact Bernoulli) ===");
kum("kappa_3^{w=1}", 3, chim3, 1, [1,2,3,4,5,6]);
kum("kappa_3^{w=2}", 3, triv,  2, [1,2,3,4,5,6]);
kum("kappa_2^{w=1}", 2, chim4, 1, [1,2,3,4,5,6,7,8]);
kum("kappa_2^{w=1,chi12}", 2, chim3, 1, [1,2,3,4,5,6,7,8]);
kum("kappa_2^{w=2}", 2, triv,  2, [1,2,3,4,5,6,7,8]);
kum("kappa_5^{w=2}", 5, chi5,  2, [1,2,3,4]);
print();
print("=== criterion FAILS: expect no p-adic convergence ===");
nocv("C   (10,3,9)",   2, rowR2(10,3,9,150));
nocv("B   (9,3,27)",   2, rowR2(9,3,27,150));
nocv("gamma Apery",    2, rowR3(17,5,1,150));
nocv("gamma Apery",    3, rowR3(17,5,1,150));
nocv("delta (7,3,81)", 2, rowR3(7,3,81,150));
nocv("alpha Domb",     3, rowR3(10,4,64,150));
nocv("eta (11,5,125)", 2, rowR3(11,5,125,150));
nocv("cusp L(f,2)",    3, rowCusp(150));
nocv("Domb* cuspidal", 2, rowDombStar(150));
nocv("Domb* cuspidal", 3, rowDombStar(150));
print();
print("=== P-strand: convergence along towers n = a p^s (criterion fails) ===");
tow("C (10,3,9)", 2, rowR2(10,3,9,600), [1,3,5]);
tow("Domb", 3, rowR3(10,4,64,600), [1,2,4]);
print();
print("TOTAL time ",gettime()," ms");
}
quit
