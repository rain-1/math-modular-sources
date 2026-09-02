/* 44_lucas.gp -- Tasks (3),(4),(5): Lucas-type congruences, Dwork/Beukers
   supercongruences, and the truncated-Dwork / P^{(p-1)/2} attack, for
   Cooper's rows s7,s10,s18 (k=1,2,3).  Exact integer arithmetic throughout. */
read("40_core.gp");

N = 600;              /* main range for A_n, L_n, a_n */
PL = [2,3,5,7,11,13,17,19,23];        /* main prime list for (3a)-(3e),(3g) */
PL43 = [2,3,5,7,11,13,17,19,23,29,31,37,41,43]; /* for (3f) */

print("################ Task (3): Lucas-type congruences, N=",N," ################");

{
for(k=1,3,
  my(nm=NAM[k]);
  print("\n========== row ",nm," (k=",k,") ==========");
  my(A = AVEC(k,N));            /* A[n+1] = A_n, n=0..N */
  my(a = AJ(k,N));               /* a[j+1] = a_j, j=0..N */
  my(Lv = vector(N, n, A[n]/n)); /* Lv[n] = L_n = A_{n-1}/n, n=1..N  (Lv index n -> n) */

  /* ---------- (3a) A_{pj+r} == A_j*A_r (mod p), 0<=r<p ---------- */
  print("-- (3a) A_{pj+r} == A_j A_r (mod p) --");
  for(pi=1,#PL, my(p=PL[pi], fail=0, tot=0, firstfail=0);
    for(j=0,N\p,
      for(r=0,min(p-1,N-p*j),
        my(n=p*j+r);
        if(n>N, next);
        tot++;
        if((A[n+1]-A[j+1]*A[r+1])%p != 0,
           fail++; if(firstfail==0, firstfail=[j,r,n]))
      )
    );
    print("  p=",p,": fails=",fail,"/",tot, if(fail>0, Str(" first fail (j,r,n)=",firstfail), ""))
  );

  /* ---------- (3b) L_{pj+r} == L_j*L_r (mod p) where defined ---------- */
  print("-- (3b) L_{pj+r} == L_j L_r (mod p), j,r>=1 --");
  for(pi=1,#PL, my(p=PL[pi], fail=0, tot=0, firstfail=0);
    for(j=1,N\p,
      for(r=1,min(p-1,N-p*j),
        my(n=p*j+r);
        if(n>N, next);
        tot++;
        if((Lv[n]-Lv[j]*Lv[r])%p != 0,
           fail++; if(firstfail==0, firstfail=[j,r,n]))
      )
    );
    print("  p=",p,": fails=",fail,"/",tot, if(fail>0, Str(" first fail (j,r,n)=",firstfail), ""))
  );

  /* ---------- (3c) a_{pj+r} == a_j*a_r (mod p) -- control, expect fail ---------- */
  print("-- (3c) a_{pj+r} == a_j a_r (mod p)  [control] --");
  for(pi=1,#PL, my(p=PL[pi], fail=0, tot=0, firstfail=0);
    for(j=0,N\p,
      for(r=0,min(p-1,N-p*j),
        my(n=p*j+r);
        if(n>N, next);
        tot++;
        if((a[n+1]-a[j+1]*a[r+1])%p != 0,
           fail++; if(firstfail==0, firstfail=[j,r,n]))
      )
    );
    print("  p=",p,": fails=",fail,"/",tot, if(fail>0, Str(" first fail (j,r,n)=",firstfail), ""))
  );

  /* ---------- (3d) Dwork/Beukers: A_{mp^s} == A_{mp^{s-1}} (mod p^{3s}), s=1,2 ---------- */
  print("-- (3d) A_{m p^s} == A_{m p^{s-1}} (mod p^{3s}), s=1,2 --");
  for(pi=1,#PL, my(p=PL[pi]);
    my(fail1=0,tot1=0,fail2=0,tot2=0,ff1=0,ff2=0);
    for(m=1,N\p,
      tot1++;
      if((A[m*p+1]-A[m+1]) % p^3 != 0, fail1++; if(ff1==0,ff1=m))
    );
    for(m=1,N\p^2,
      tot2++;
      if((A[m*p^2+1]-A[m*p+1]) % p^6 != 0, fail2++; if(ff2==0,ff2=m))
    );
    print("  p=",p,": s=1 fails=",fail1,"/",tot1,(if(fail1>0, Str(" firstfail m=",ff1), "")),
          "   s=2 fails=",fail2,"/",tot2,(if(fail2>0, Str(" firstfail m=",ff2), "")))
  );

  /* ---------- (3e) largest e with A_{mp}==A_m (mod p^e) for all m<=40 ---------- */
  print("-- (3e) largest e: A_{mp} == A_m (mod p^e) for all m<=40 (m capped by N) --");
  for(pi=1,#PL, my(p=PL[pi], mmax=min(40,N\p));
    my(e=0);
    while(1,
      my(ok=1);
      for(m=1,mmax, if((A[m*p+1]-A[m+1]) % p^(e+1) != 0, ok=0; break));
      if(ok, e++, break);
      if(e>30, break) /* safety cap */
    );
    print("  p=",p,": mmax used=",mmax,"  largest e = ",e)
  );

  /* ---------- (3f) p | A_{p-1}? and v_p(A_{p-1}), p<=43 ---------- */
  print("-- (3f) v_p(A_{p-1}), p<=43 --");
  my(Abig = if(N>=42, A, AVEC(k,42)));
  for(pi=1,#PL43, my(p=PL43[pi]);
    if(p-1<=N, my(v=valuation(A[p], p)); print("  p=",p,": A_{p-1}=A_",p-1,"  v_p=",v),
       my(A2=AVEC(k,p-1)); my(v=valuation(A2[p],p)); print("  p=",p,": A_{p-1}=A_",p-1," v_p=",v))
  );

  /* ---------- (3g) truncated-Dwork == test (3a); report first failing (p,j,r) ---------- */
  print("-- (3g) truncated Dwork [same data as (3a)]: first failing (p,j,r) per row --");
  my(anyfail=0);
  for(pi=1,#PL, my(p=PL[pi]);
    for(j=0,N\p,
      for(r=0,min(p-1,N-p*j),
        my(n=p*j+r);
        if(n>N, next);
        if((A[n+1]-A[j+1]*A[r+1])%p != 0,
           print("  FIRST FAIL: p=",p," j=",j," r=",r," n=",n); anyfail=1; break(3))
      )
    )
  );
  if(!anyfail, print("  no failures found for any p in PL over n<=",N));
)
};

print("\n\n################ Task (4): b_j := [x^j] l(x)/(x F(x)) ################");
{
for(k=1,3,
  my(nm=NAM[k]);
  print("\n========== row ",nm," (k=",k,") ==========");
  my(BB = LSER(k,N)/('x*FSER(k,N)));
  my(bvec = vector(N+1,j, polcoeff(BB,j-1))); /* bvec[j+1]=b_j, j=0..N */
  print("-- (4a) b_{pj+r} == b_j b_r (mod p) --");
  for(pi=1,#PL, my(p=PL[pi], fail=0, tot=0, firstfail=0);
    for(j=0,N\p,
      for(r=0,min(p-1,N-p*j),
        my(n=p*j+r);
        if(n>N, next);
        tot++;
        if((bvec[n+1]-bvec[j+1]*bvec[r+1])%p != 0,
           fail++; if(firstfail==0, firstfail=[j,r,n]))
      )
    );
    print("  p=",p,": fails=",fail,"/",tot, if(fail>0, Str(" first fail (j,r,n)=",firstfail), ""))
  );
  print("-- (4b) b_{pj+p-1} == psi(p) b_j (mod p)  [control, expect FAIL] --");
  for(pi=1,#PL, my(p=PL[pi], fail=0, tot=0, firstfail=0, ps=psival(k,p));
    for(j=0,N\p,
      my(n=p*j+p-1);
      if(n>N, next);
      tot++;
      if((bvec[n+1]-ps*bvec[j+1])%p != 0,
         fail++; if(firstfail==0, firstfail=[j,n]))
    );
    print("  p=",p,": fails=",fail,"/",tot, if(fail>0, Str(" first fail (j,n)=",firstfail), ""))
  );
)
};

print("\n\n################ Task (5): P^{(p-1)/2} attack ################");
PL31 = [2,3,5,7,11,13,17,19,23,29,31];
{
for(k=1,3,
  my(nm=NAM[k]);
  print("\n========== row ",nm," (k=",k,") ==========");
  my(BB = LSER(k,N)/('x*FSER(k,N)));
  my(bvec = vector(N+1,j, polcoeff(BB,j-1)));
  my(P = PSER(k,N));
  print("-- (5a) [x^{pj+p-1}](b*P^{(p-1)/2}) == psi(p) b_j  (mod p) --");
  for(pi=1,#PL31, my(p=PL31[pi]);
    my(e=(p-1)/2);
    my(Pe = (P^e + O('x^(N+1))));
    my(prod = BB*Pe + O('x^(N+1)));
    my(cvec = vector(N+1,j,polcoeff(prod,j-1)));
    my(ps=psival(k,p), fail=0, tot=0, firstfail=0);
    for(j=0,N\p,
      my(n=p*j+p-1);
      if(n>N, next);
      tot++;
      if((cvec[n+1]-ps*bvec[j+1])%p != 0,
         fail++; if(firstfail==0, firstfail=[j,n]))
    );
    print("  p=",p,": fails=",fail,"/",tot, if(fail>0, Str(" first fail (j,n)=",firstfail), ""));
  );
  print("-- (5b) same test mod p^2 --");
  for(pi=1,#PL31, my(p=PL31[pi]);
    my(e=(p-1)/2);
    my(Pe = (P^e + O('x^(N+1))));
    my(prod = BB*Pe + O('x^(N+1)));
    my(cvec = vector(N+1,j,polcoeff(prod,j-1)));
    my(ps=psival(k,p), fail=0, tot=0, firstfail=0);
    for(j=0,N\p,
      my(n=p*j+p-1);
      if(n>N, next);
      tot++;
      if((cvec[n+1]-ps*bvec[j+1])%p^2 != 0,
         fail++; if(firstfail==0, firstfail=[j,n]))
    );
    print("  p=",p,": fails=",fail,"/",tot, if(fail>0, Str(" first fail (j,n)=",firstfail), ""));
  );
)
};

print("\nDONE.");
quit;
