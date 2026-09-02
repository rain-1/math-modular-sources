\\ 37_level28.gp -- Task (4): weakly holomorphic weight-5/2 forms at level 28.
\\ Space W_r = { g / Delta(4 tau)^r : g in M_{5/2+12r}(Gamma_0(28)) }.
\\ Look for f in the Kohnen plus space with principal part exactly q^{-3} and a(0)=0,
\\ then test  a(m^2) = m * beta_{s7}(m)  (Cooper row s7, psi=1).
default(parisize, 12000000000);
default(timer,1);
NC = 420;                       \\ number of q-coefficients kept
BS7 = read("beta_s7.txt");      \\ beta_{s7}(1..) from 01_data.gp
print("beta_s7(1..20) = ", vector(20,i,BS7[i]));
D4 = 'q^4;
{ D4 = 'q^4*prod(n=1, (NC+8)\4, (1-'q^(4*n))^24) + O('q^(NC+10)); }
{ run(r) =
  my(w, MM, B, d, S, L, A, e, sol, ker, NCC);
  w = 5/2 + 12*r;
  print("");
  print("############ r = ", r, " ,  weight ", w, " ############");
  MM = mfinit([28,w]);
  d = mfdim(MM);
  print("  dim M_",w,"(Gamma_0(28)) = ", d);
  B = mfbasis(MM);
  \\ S[i] = B[i]/Delta(4tau)^r  as a coefficient vector indexed n = -4r .. NC
  S = vector(d);
  for(i=1,d,
    my(ser = Ser(mfcoefs(B[i], NC+4*r+4), 'q) + O('q^(NC+4*r+5)));
    ser = ser/D4^r;
    S[i] = vector(NC+4*r+1, t, polcoeff(ser, t-1-4*r)));
  \\ index helper: S[i][n + 4r + 1] = coefficient of q^n
  \\ conditions
  L = [];
  for(n=-4*r, -1, if(n != -3, L = concat(L,[n])));
  L = concat(L, [0]);
  for(n=1, NC, if(n%4==2 || n%4==3, L = concat(L,[n])));
  print("  #linear conditions = ", #L, "   #unknowns = ", d);
  A = matrix(#L+1, d, ii, jj, if(ii<=#L, S[jj][L[ii]+4*r+1], S[jj][-3+4*r+1]));
  e = vectorv(#L+1); e[#L+1] = 1;
  sol = 0;
  iferr(sol = matinverseimage(A, e), err, print("  matinverseimage error: ", err));
  if(type(sol)!="t_COL" || #sol==0,
     print("  ==> NO element of W_",r," has principal part exactly q^-3 in the plus space."),
     ker = matker(A);
     print("  ==> solution EXISTS; solution space dimension = ", #ker);
     my(av = vector(NC+4*r+1));
     for(t=1,NC+4*r+1, av[t] = sum(i=1,d, sol[i]*S[i][t]));
     my(ac(n) = av[n+4*r+1]);
     print("  principal part: ", vector(4*r, t, Str("a(",-4*r+t-1,")=",ac(-4*r+t-1))), "  a(0)=", ac(0));
     print("  a(n), n=1..40 : ", vector(40,n,ac(n)));
     print("  --- TEST  a(m^2) = m * beta_s7(m) ---");
     my(nagree=0, ndis=0);
     for(m=1, 20, my(x=ac(m^2), y=m*BS7[m]);
       print("    m=",m,"  a(m^2)=",x,"   m*beta_s7=",y,"   ", if(x==y,"MATCH","differ"));
       if(x==y, nagree++, ndis++));
     print("  matches: ", nagree, " / 20");
     if(#ker>0,
       print("  --- with the ", #ker, "-dim ambiguity: can it be fixed to match? ---");
       my(KA = matrix(20, #ker, ii, jj, sum(i=1,d, ker[i,jj]*S[i][ii^2+4*r+1])),
          tgt = vectorv(20, ii, ii*BS7[ii] - ac(ii^2)), fx);
       fx = 0;
       iferr(fx = matinverseimage(KA, tgt), err2, print("    (no solution: ", err2, ")"));
       if(type(fx)=="t_COL" && #fx>0,
          print("    YES: the ambiguity can be tuned; correction = ", fx),
          print("    NO: no choice in the ", #ker, "-dim ambiguity reproduces beta_s7 on m=1..20")));
  );
}
run(1);
run(2);
quit;
