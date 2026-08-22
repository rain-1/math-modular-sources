/* towers2.gp -- tower Frobenius eigenvalues for the NON-MODULAR rows:
   Zudilin's Catalan row and the Brown-Zudilin zeta(5) cellular row.
   Measures  rhoA_s = A(a p^{s+1})/A(a p^s)   (should -> 1, Dwork)
             rhoB_s = B(a p^{s+1})/B(a p^s)   (should -> eps * p^{-(w+1)})
   by reporting v_p(rhoA_s - 1) and v_p(u * p^{w+1} * rhoB_s - 1) for u = +-1. */
default(parisizemax, 6000000000);
read("lattice/euler_criterion/rows.gp");

c0(n) = 2*(2*n+1)*(41218*n^3-48459*n^2+20010*n-2871)*(n+1)^5;
c1(n) = -(97604224*n^9 + 178061760*n^8 + 72005308*n^7 - 48634688*n^6 - 39076836*n^5 + 2622730*n^4 + 7581006*n^3 + 920112*n^2 - 543402*n - 120582);
c2(n) = -2*n*(3874492*n^8 - 2617900*n^7 - 3144314*n^6 + 2947148*n^5 + 647130*n^4 - 1182926*n^3 + 115771*n^2 + 170716*n - 44541);
c3(n) = n*(41218*n^3+75195*n^2+46746*n+9898)*(n-1)^5;
{ bzrun(v0,v1,v2,M) = my(v=vector(M+1)); v[1]=v0; v[2]=v1; v[3]=v2;
  for(n=2,M-1, v[n+2] = -(c1(n)*v[n+1] + c2(n)*v[n] + c3(n)*v[n-1])/c0(n)); v; }

vv(x,p) = if(x==0, "inf", valuation(x,p));
{ tower(nm, A, B, M, p, a, w1) =
  my(smax = floor(log(M/a)/log(p)), sA=[], sB=[], sBm=[]);
  for(s=0, smax-1,
    my(n0=a*p^s, n1=a*p^(s+1));
    if(n1>M, break);
    my(rA = A[n1+1]/A[n0+1], rB = if(B[n0+1]==0, 0, B[n1+1]/B[n0+1]));
    sA = concat(sA, [vv(rA-1,p)]);
    if(rB!=0,
      sB  = concat(sB,  [vv(p^w1*rB-1,p)]);
      sBm = concat(sBm, [vv(-p^w1*rB-1,p)])));
  print(nm, "  p=",p," a=",a,"  v(rhoA-1)=",sA);
  print("        v(+p^",w1,"rhoB-1)=",sB,"   v(-p^",w1,"rhoB-1)=",sBm); }

MZ = 2200; Z = rowZud(MZ);
print("=== Zudilin Catalan row, w+1 = 2 ===");
{ for(j=1,4, my(p=[3,5,7,11][j]); for(a=1,2, tower("Zud", Z[1], Z[2], MZ, p, a, 2))); }

MB = 900; QB = bzrun(1,21,2989,MB); PB = bzrun(0, 87/4, 1190161/384, MB);
print("=== Brown-Zudilin zeta(5) row, w+1 = 5 ===");
{ for(j=1,4, my(p=[2,3,5,7][j]); for(a=1,2, tower("BZ ", QB, PB, MB, p, a, 5))); }
print("=== BZ zeta(3) direction (Phat), w+1 = 3 ? ===");
PH = bzrun(0, 101/4, 344923/96, MB);
{ for(j=1,2, my(p=[2,3][j]); for(a=1,2, tower("BZh", QB, PH, MB, p, a, 3))); }
