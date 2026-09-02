\\ 38_wh.gp -- Task (4).  V_r(N) := { f = g/Delta(4tau)^r : g in M_{5/2+12r}(Gamma_0(4N)),
\\   f in the Kohnen plus space, a(n)=0 for n<0 with n != -3, a(0)=0 }.
\\ No scale is imposed (at level 1 the answer is f4a = (1/64)q^-3 + q - 506 q^4 + ...).
\\ Question: does V_r contain a form with a(m^2) = m*beta(m) for the target row?
\\ Output goes line by line to 38_wh_out.txt (write() flushes immediately).
default(parisize, 14000000000);
OUT = "38_wh_out.txt";
say(s) = write(OUT, s);
NC = 410;
BS7 = read("beta_s7.txt");
{ NNc = 60;
  E4s = 1 + 240*sum(n=1, NNc-1, sigma(n,3)*'q^n) + O('q^NNc);
  E6s = 1 - 504*sum(n=1, NNc-1, sigma(n,5)*'q^n) + O('q^NNc);
  DELs = (E4s^3 - E6s^2)/1728;
  Aa = vector(50, m, polcoeff(DELs/E4s^2, m));
  CPa = vector(50, n, Aa[n]/n);
  BCTL = vector(50);
  for(n=1,50, my(s=0); fordiv(n,d, s+=moebius(d)*CPa[n/d]); BCTL[n]=s);
}
{ run(N, r, TGT, nam) =
  my(w, MM, B, d, S, L, A, ker, PRE, NT, dv, VB, KA, tgt, sol, kmax, D4);
  w = 5/2 + 12*r;
  PRE = NC + 4*r + 30;
  say(Str("######## ", nam, " : level ", 4*N, ", r=", r, ", weight ", w, " ########"));
  MM = mfinit([4*N, w]);
  d = mfdim(MM);
  say(Str("  dim M_",w,"(Gamma_0(",4*N,")) = ", d));
  if(d==0, return(0));
  B = mfbasis(MM);
  D4 = 'q^4*prod(n=1, (PRE+8)\4, (1-'q^(4*n))^24) + O('q^(PRE+6));
  S = vector(d);
  for(i=1,d, my(ser = Ser(mfcoefs(B[i], PRE), 'q) + O('q^(PRE+1)));
    ser = ser/D4^r;
    S[i] = vector(NC+4*r+1, t, polcoeff(ser, t-1-4*r)));
  say("  q-expansions built");
  L = [];
  for(n=-4*r, -1, if(n != -3, L = concat(L,[n])));
  L = concat(L, [0]);
  for(n=1, NC, if(n%4==2 || n%4==3, L = concat(L,[n])));
  A = matrix(#L, d, ii, jj, S[jj][L[ii]+4*r+1]);
  ker = matker(A);
  dv = #ker;
  say(Str("  #homogeneous conditions = ", #L, "   dim V_r = ", dv));
  if(dv==0, say("  ==> V_r = 0, nothing to test"); return(0));
  VB = vector(dv, jj, vector(NC+4*r+1, t, sum(i=1,d, ker[i,jj]*S[i][t])));
  say(Str("  pole orders of the V_r basis: ", vector(dv, jj, my(v=-4*r); while(v<=NC && VB[jj][v+4*r+1]==0, v++); v)));
  say(Str("  a(-3) on the V_r basis: ", vector(dv, jj, VB[jj][-3+4*r+1])));
  NT = min(20, sqrtint(NC));
  KA  = matrix(NT, dv, ii, jj, VB[jj][ii^2+4*r+1]);
  tgt = vectorv(NT, ii, ii*TGT[ii]);
  sol = 0;
  iferr(sol = matinverseimage(KA, tgt), er, sol = 0);
  if(type(sol)=="t_COL",
     say(Str("  ==> *** EXACT MATCH on m=1..", NT, " *** coeff vector = ", sol));
     my(fin = vector(NC+4*r+1, t, sum(jj=1,dv, sol[jj]*VB[jj][t])));
     say(Str("      a(n), n=", -4*r, "..40 : ", vector(4*r+41, t, fin[t]))),
     say(Str("  ==> NO f in V_r has a(m^2)=m*beta(m) for m=1..", NT, "  (rank ", matrank(KA), " of the ", NT, "x", dv, " system)"));
     kmax = 0;
     for(K=1,NT, my(KB=matrix(K,dv,ii,jj,KA[ii,jj]), tb=vectorv(K,ii,tgt[ii]), z=0);
       iferr(z=matinverseimage(KB,tb), e3, z=0);
       if(type(z)=="t_COL", kmax=K, break));
     say(Str("      longest matchable prefix m=1..k : k = ", kmax)));
  1;
}
say("############ CONTROL: level 4 (the answer must be proportional to f4a) ############");
run(1, 1, BCTL, "CONTROL level 4 r=1");
say("");
say("############ TARGET: level 28, Cooper row s7 (psi=1) ############");
run(7, 1, BS7, "s7 r=1");
run(7, 2, BS7, "s7 r=2");
say("DONE");
quit;
