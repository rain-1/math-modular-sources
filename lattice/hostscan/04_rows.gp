/* 04_rows.gp -- build the row on each unit-lambda_2 Fricke host and compute
   its arithmetic: 3-term recurrence, free integration, denominator exponent k,
   Apery limit xi (high precision).   Theorem 3.4 shape.  */
default(parisizemax, 8000000000);
default(realprecision, 120);
NQ = 140;          /* q-series order */
NA = 110;          /* number of a_n peeled */
NB = 400;          /* recurrence iterations for xi and for k */

/* eta product series: prod_d (q^{d/24} prod(1-q^{dn}))^{r_d}, with q^{sum d r_d/24}=q */
{useries(nn, dv, r, nq) =
  my(s = 1 + O(q^nq));
  for(t=1,#dv, my(d=dv[t], e=r[t]);
    if(e!=0, s *= eta(q^d + O(q^nq))^e));
  q*s;}

/* F = D log u = 1 - sum_d r_d d sum_{n>=1} sigma_1(n) q^{dn} */
{Fseries(nn, dv, r, nq) =
  my(s = 1 + O(q^nq));
  for(t=1,#dv, my(d=dv[t], e=r[t]);
    if(e!=0, s -= e*d*sum(n=1, (nq-1)\d, sigma(n)*q^(d*n)) ));
  s + O(q^nq);}

/* peel F = sum a_n x^n */
{peel(Fs, xs, na, nq) =
  my(a = vector(na+1), G = Fs, xp = 1 + O(q^nq));
  for(n=0, na,
    my(c = polcoeff(G, n));
    a[n+1] = c;
    G = G - c*xp2(xs, n, nq);
  );
  a;}
{xp2(xs, n, nq) = if(n==0, 1+O(q^nq), xs^n);}

/* better peel: iterative */
{peel2(Fs, xs, na, nq) =
  my(a = vector(na+1), G = Fs, xp = 1 + O(q^nq));
  for(n=0, na,
    my(c = polcoeff(G, n));
    a[n+1] = c;
    if(c != 0, G = G - c*xp);
    xp = xp*xs;
  );
  a;}

/* fit sum_{j=0}^{2} p_j(n) a_{n+j} = 0, deg p_j <= dg ; return kernel basis */
{fitrec(av, dg) =
  my(nv = 3*(dg+1), rows = List());
  for(n=0, #av-1-2-1,
    my(row = vector(nv));
    for(j=0,2, for(e=0,dg, row[j*(dg+1)+e+1] = n^e * av[n+j+1]));
    listput(rows, row));
  my(mm = matconcat(Vec(rows)~));
  matker(mm);}

{run(nn, dv, r, C, B) =
  my(us = useries(nn,dv,r,NQ), Fs = Fseries(nn,dv,r,NQ));
  my(xs = us/(1 + B*us + C*us^2));
  my(a = peel2(Fs, xs, NA, NQ));
  /* check integrality */
  my(intok = 1); for(t=1,#a, if(denominator(a[t])!=1, intok=0; break));
  my(ker = fitrec(a, 3));
  my(rr = matsize(ker)[2]);
  print("N=",nn," C=",C," B=",B," int=",intok," recdim=",rr," a=",vector(min(8,#a),i,a[i]));
  if(rr != 1, print("   ** recurrence dim ",rr," -- not a clean 3-term/deg3 row"); return(0));
  my(kv = ker[,1]);
  /* normalise */
  my(P = vector(3, j, sum(e=0,3, kv[(j-1)*4+e+1]*x^e)));
  my(g = content(concat(Vec(P[1]),concat(Vec(P[2]),Vec(P[3]))))); P = vector(3,j,P[j]/g);
  print("   P0(n)=",P[1],"  P1(n)=",P[2],"  P2(n)=",P[3]);
  /* char roots from leading coefficients */
  my(cp = sum(j=1,3, polcoeff(P[j],3)*y^(j-1)));
  print("   charpoly=",cp,"  roots=",polroots(cp));
  /* free integration test */
  my(freeint = 1); for(n=1,60, if(a[n+1] % (n+1) != 0, freeint = 0; break));
  print("   free integration (n+1)|A_n: ", if(freeint,"YES","no"));
  /* companion */
  my(bb = vector(NB+2)); bb[1]=0; bb[2]=1;
  my(aa = vector(NB+2)); for(t=1,min(#a,NB+2), aa[t]=a[t]);
  /* extend a by recurrence too, and compute b */
  for(n=1, NB,
    my(p0 = subst(P[1],x,n), p1 = subst(P[2],x,n), p2 = subst(P[3],x,n));
    if(p2==0, print("   ** leading coeff vanishes at n=",n); break);
    if(n+2 <= NB+2,
      aa[n+2] = -(p0*aa[n] + p1*aa[n+1])/p2;
      bb[n+2] = -(p0*bb[n] + p1*bb[n+1])/p2));
  /* check aa matches a */
  my(amatch=1); for(t=1,min(#a,NB+2), if(aa[t]!=a[t], amatch=0; break));
  print("   recurrence reproduces a_n: ", if(amatch,"YES","NO"));
  /* denominator exponent k */
  my(dn = 1, kmax = 0, kk);
  for(n=1, 150,
    dn = lcm(dn, n);
    my(bn = bb[n+1]);
    if(bn != 0, my(de = denominator(bn), kj = 0, t = de);
      while(t > 1, kj++; t = t/gcd(t, dn); if(kj>9, break));
      if(kj > kmax, kmax = kj)));
  print("   denominator exponent k = ", kmax);
  /* Apery limit */
  my(xi = bb[NB+1]/aa[NB+1]*1.0, xi2 = bb[NB]/aa[NB]*1.0);
  print("   xi ~ ", xi, "   (prev ", xi2, ")");
  [P, kmax, freeint, xi];
}
