\\ 39_v28.gp -- level 28, weight 5/2 weakly holomorphic, TWO support conditions.
\\   V_r  : scalar Kohnen plus,  a(d) = 0 unless d = 0,1 mod 4
\\   V'_r : the s7 support law,  a(d) = 0 unless -3d is a square mod 28
\\          (equivalently d = 0,1 mod 4 AND d = 0,1,2,4 mod 7) -- this is the support
\\          observed for the twisted CM traces of the s7 row.
\\ In each, impose the only pole at oo to be q^{-3} and a(0)=0, then test
\\   a(m^2) = m * beta_s7(m).
\\ Output to 39_v28_out.txt.
default(parisize, 12000000000);
OUT = "39_v28_out.txt";
say(s) = write(OUT, s);
NC = 1700;
BS7 = read("beta_s7.txt");
{ supp(d, mode) =
  \\ returns 1 if the coefficient a(d) is ALLOWED, 0 if it must vanish
  if(mode==0, return(d%4==0 || d%4==1));
  if((d%4!=0 && d%4!=1), return(0));
  my(s = ((-3*d)%7+7)%7);
  (s==0 || s==1 || s==2 || s==4);
}
{ run(N, r, mode, TGT, nam) =
  my(w, MM, B, d, S, L, A, ker, PRE, NT, dv, VB, KA, tgt, sol, kmax, D4, tau);
  w = 5/2 + 12*r;
  PRE = NC + 4*r + 30;
  say(Str("######## ", nam, " : level ", 4*N, ", r=", r, ", weight ", w,
          ", support mode ", mode, if(mode==0, " (plus only)", " (plus AND -3d square mod 28)"), " ########"));
  MM = mfinit([4*N, w]);
  d = mfdim(MM);
  say(Str("  dim M_",w,"(Gamma_0(",4*N,")) = ", d));
  if(d==0, return(0));
  B = mfbasis(MM);
  tau = mfcoefs(mfDelta(), (PRE+8)\4);
  D4 = Ser(vector(PRE+6, t, if((t-1)%4==0 && (t-1)\4 <= #tau-1, tau[(t-1)\4+1], 0)), 'q) + O('q^(PRE+6));
  S = vector(d);
  for(i=1,d, my(ser = Ser(mfcoefs(B[i], PRE), 'q) + O('q^(PRE+1)));
    ser = ser/D4^r;
    S[i] = vector(NC+4*r+1, t, polcoeff(ser, t-1-4*r)));
  say("  q-expansions built");
  L = [];
  for(n=-4*r, -1, if(n != -3, L = concat(L,[n])));
  L = concat(L, [0]);
  for(n=1, NC, if(!supp(n, mode), L = concat(L,[n])));
  A = matrix(#L, d, ii, jj, S[jj][L[ii]+4*r+1]);
  ker = matker(A);
  dv = #ker;
  say(Str("  #homogeneous conditions = ", #L, "   dim = ", dv));
  if(dv==0, say("  ==> the space is 0"); return(0));
  VB = vector(dv, jj, vector(NC+4*r+1, t, sum(i=1,d, ker[i,jj]*S[i][t])));
  say(Str("  pole orders: ", vector(dv, jj, my(v=-4*r); while(v<=NC && VB[jj][v+4*r+1]==0, v++); v)));
  NT = min(40, sqrtint(NC));
  KA  = matrix(NT, dv, ii, jj, VB[jj][ii^2+4*r+1]);
  tgt = vectorv(NT, ii, ii*TGT[ii]);
  sol = 0;
  iferr(sol = matinverseimage(KA, tgt), er, sol = 0);
  if(type(sol)=="t_COL" && #sol>0,
     say(Str("  ==> *** EXACT MATCH on m=1..", NT, " ***"));
     my(fin = vector(NC+4*r+1, t, sum(jj=1,dv, sol[jj]*VB[jj][t])));
     say(Str("      a(n), n=", -4*r, "..120 : ", vector(4*r+121, t, fin[t]))),
     say(Str("  ==> NO element has a(m^2)=m*beta_s7(m) for m=1..", NT,
             "  (rank ", matrank(KA), " of the ", NT, "x", dv, " system)"));
     kmax = 0;
     for(K=1,NT, my(KB=matrix(K,dv,ii,jj,KA[ii,jj]), tb=vectorv(K,ii,tgt[ii]), z=0);
       iferr(z=matinverseimage(KB,tb), e3, z=0);
       if(type(z)=="t_COL" && #z>0, kmax=K, break));
     say(Str("      longest matchable prefix m=1..k : k = ", kmax,
             "   (k = dim means the system behaves generically, i.e. no hidden solution)")));
  1;
}
say("############ level 28, s7 : plus-only vs plus+mod-7 support ############");
run(7, 1, 1, BS7, "s7 r=1, s7-support");
run(7, 2, 1, BS7, "s7 r=2, s7-support");
say("DONE");
quit;
