\\ 38_wh.gp -- Task (4), corrected & with a level-4 CONTROL.
\\ W_r(N) = { g/Delta(4tau)^r : g in M_{5/2+12r}(Gamma_0(4N)) }, weight 5/2, poles at cusps.
\\ Inside W_r(N) impose: Kohnen plus (a(n)=0 for n=2,3 mod 4), a(n)=0 for n<0, n!=-3,
\\ a(-3)=1, a(0)=0.  Then compare a(m^2) with m*beta(m) for the target row.
default(parisize, 14000000000);
default(timer,1);
NC   = 420;
BS7  = read("beta_s7.txt");
\\ control target: beta for F4a = Delta/E4^2 (level 4), where the answer is 64 f4a
{ NNc = 60;
  E4s = 1 + 240*sum(n=1, NNc-1, sigma(n,3)*'q^n) + O('q^NNc);
  E6s = 1 - 504*sum(n=1, NNc-1, sigma(n,5)*'q^n) + O('q^NNc);
  DELs = (E4s^3 - E6s^2)/1728;
  Aa = vector(50, m, polcoeff(DELs/E4s^2, m));
  CPa = vector(50, n, Aa[n]/n);
  BCTL = vector(50); for(n=1,50, my(s=0); fordiv(n,d, s+=moebius(d)*CPa[n/d]); BCTL[n]=s);
}
{ run(N, r, TGT, tgtscale, nam) =
  my(w, MM, B, d, S, L, A, e, sol, ker, PRE, ac);
  w = 5/2 + 12*r;
  PRE = NC + 4*r + 40;
  print("");
  print("######## ", nam, " : level 4*", N, "=", 4*N, ", r=", r, ", weight ", w, " ########");
  MM = mfinit([4*N, w]);
  d = mfdim(MM);
  print("  dim M_",w,"(Gamma_0(",4*N,")) = ", d);
  if(d==0, print("  empty"); return(0));
  B = mfbasis(MM);
  my(D4 = 'q^4*prod(n=1, (PRE+8)\4, (1-'q^(4*n))^24) + O('q^(PRE+6)));
  S = vector(d);
  for(i=1,d,
    my(ser = Ser(mfcoefs(B[i], PRE), 'q) + O('q^(PRE+1)));
    ser = ser/D4^r;
    S[i] = vector(NC+4*r+1, t, polcoeff(ser, t-1-4*r)));
  L = [];
  for(n=-4*r, -1, if(n != -3, L = concat(L,[n])));
  L = concat(L, [0]);
  for(n=1, NC, if(n%4==2 || n%4==3, L = concat(L,[n])));
  print("  #conditions = ", #L+1, "   #unknowns = ", d);
  A = matrix(#L+1, d, ii, jj, if(ii<=#L, S[jj][L[ii]+4*r+1], S[jj][-3+4*r+1]));
  e = vectorv(#L+1); e[#L+1] = 1;
  sol = 0;
  iferr(sol = matinverseimage(A, e), err, print("  no solution (matinverseimage): ", err); sol=0);
  if(type(sol)!="t_COL",
     print("  ==> NO element with principal part exactly q^-3 in the plus space of W_",r,"."); return(0));
  ker = matker(A);
  print("  ==> EXISTS.  ambiguity dim = ", #ker, "  (= dim of plus forms in W_r holomorphic at oo, vanishing there)");
  my(av = vector(NC+4*r+1, t, sum(i=1,d, sol[i]*S[i][t])));
  ac(n) = av[n+4*r+1];
  print("  a(n) n=-4r..6 : ", vector(4*r+7, t, Str(t-1-4*r,":",ac(t-1-4*r))));
  my(NT = min(20, sqrtint(NC)));
  print("  --- TEST a(m^2) =? ", tgtscale, "*m*beta(m), m=1..", NT, " ---");
  my(nag=0);
  for(m=1,NT, if(ac(m^2)==tgtscale*m*TGT[m], nag++));
  print("    exact matches for the particular solution: ", nag, "/", NT);
  if(#ker>0,
    my(KA = matrix(NT, #ker, ii, jj, sum(i=1,d, ker[i,jj]*S[i][ii^2+4*r+1])),
       tgt = vectorv(NT, ii, tgtscale*ii*TGT[ii] - ac(ii^2)), fx=0);
    iferr(fx = matinverseimage(KA, tgt), er2, fx=0);
    if(type(fx)=="t_COL",
       print("    ==> AFTER tuning the ", #ker, "-dim ambiguity: EXACT MATCH on m=1..", NT);
       my(fin = vector(NC+4*r+1, t, av[t] + sum(jj=1,#ker, fx[jj]*sum(i=1,d, ker[i,jj]*S[i][t]))));
       print("    resulting a(n), n=-4r..30 : ", vector(4*r+31, t, fin[t]));
       print("    check principal part again: ", vector(4*r+1, t, Str(t-1-4*r,":",fin[t]))),
       print("    ==> NO tuning of the ", #ker, "-dim ambiguity matches on m=1..", NT,
             "  (rank ", matrank(KA), " of a ", NT, "x", #ker, " system)");
       \\ how many leading constraints CAN be matched?
       my(kmax=0);
       for(K=1,NT, my(KB=matrix(K,#ker,ii,jj,KA[ii,jj]), tb=vectorv(K,ii,tgt[ii]), z=0);
         iferr(z=matinverseimage(KB,tb), e3, z=0);
         if(type(z)=="t_COL", kmax=K, break));
       print("    largest prefix m=1..k that CAN be matched: k = ", kmax)));
  1;
}
print("################ CONTROL: level 4, target = 64*f4a (beta of Delta/E4^2) ################");
run(1, 1, BCTL, 1, "CONTROL level 4");
print("");
print("################ TARGET: level 28, Cooper row s7 ################");
run(7, 1, BS7, 1, "s7 r=1");
run(7, 2, BS7, 1, "s7 r=2");
quit;
