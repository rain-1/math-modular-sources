\\ 39_v28.gp -- level 28, weight 5/2 weakly holomorphic, DECISIVE over-determined test.
\\
\\ V'_r := { f = g/Delta(4tau)^r : g in M_{5/2+12r}(Gamma_0(28)),
\\           a(n) = 0 for every n < 0 except n = -3, a(0) = 0,
\\           a(d) = 0 for every 1 <= d <= NC that is NOT admissible },
\\ where "admissible" = ( -3d is a square mod 28 ) = ( d = 0,1 mod 4 and d = 0,1,2,4 mod 7 ).
\\ This is the support law of the twisted CM traces of Cooper's row s7 (companion strand).
\\
\\ TEST A: does some f in V'_r satisfy a(d) = kappa * c(d) for all admissible d <= NC,
\\         for a single constant kappa != 0?   (c(d) = the twisted trace table, 70_cd_s7.txt;
\\         d = 49 is degenerate there and is left unconstrained.)
\\ TEST B: does some f in V'_r satisfy the square-index dictionary a(m^2) = m*beta_s7(m)?
\\ For reference the same two tests are run in the bigger space V_r (plus condition only).
default(parisize, 12000000000);
OUT = "39_v28_out.txt";
say(s) = write(OUT, s);
NC = 400;
BS7 = read("beta_s7.txt");
CDL = read("70_cd_s7.txt");
CD = vector(NC);
DEGEN = [];
{ for(i=1, #CDL, my(e = CDL[i]);
    if(e[1] <= NC,
      if(type(e[2])=="t_STR", DEGEN = concat(DEGEN,[e[1]]), CD[e[1]] = e[2]))); }
{ adm(d) = if(d<=0, 0, if(d%4!=0 && d%4!=1, 0, my(s=((-3*d)%7+7)%7); s==0||s==1||s==2||s==4)); }
ADM = select(d->adm(d), vector(NC,i,i));
say(Str("admissible d <= ", NC, " : ", #ADM, " of them; degenerate (unconstrained): ", DEGEN));
say(Str("c(d) known for ", sum(i=1,NC, CD[i]!=0), " nonzero d"));
{ run(r, mode) =
  my(w, MM, B, d, S, L, A, ker, PRE, dv, VB, D4, tau, nam, KA, tgt, sol, kmax, cols, rows);
  w = 5/2 + 12*r;
  PRE = NC + 4*r + 30;
  nam = if(mode==0, "V_r (Kohnen plus only)", "V'_r (plus AND -3d square mod 28)");
  say("");
  say(Str("######## level 28, r=", r, ", weight ", w, " : ", nam, " ########"));
  MM = mfinit([28, w]);
  d = mfdim(MM);
  B = mfbasis(MM);
  tau = mfcoefs(mfDelta(), (PRE+8)\4);
  D4 = Ser(vector(PRE+6, t, if((t-1)%4==0 && (t-1)\4 <= #tau-1, tau[(t-1)\4+1], 0)), 'q) + O('q^(PRE+6));
  S = vector(d);
  for(i=1,d, my(ser = Ser(mfcoefs(B[i], PRE), 'q) + O('q^(PRE+1)));
    ser = ser/D4^r;
    S[i] = vector(NC+4*r+1, t, polcoeff(ser, t-1-4*r)));
  say(Str("  dim M = ", d, " ; q-expansions built"));
  L = [];
  for(n=-4*r, -1, if(n != -3, L = concat(L,[n])));
  L = concat(L, [0]);
  for(n=1, NC, if(mode==0, if(n%4==2||n%4==3, L=concat(L,[n])), if(!adm(n), L=concat(L,[n]))));
  A = matrix(#L, d, ii, jj, S[jj][L[ii]+4*r+1]);
  ker = matker(A);
  dv = #ker;
  say(Str("  #support conditions = ", #L, "   dim = ", dv));
  if(dv==0, say("  ==> space is 0"); return(0));
  VB = vector(dv, jj, vector(NC+4*r+1, t, sum(i=1,d, ker[i,jj]*S[i][t])));
  say(Str("  pole orders: ", vector(dv, jj, my(v=-4*r); while(v<=NC && VB[jj][v+4*r+1]==0, v++); v)));
  \\ ---- TEST A : a(d) = kappa c(d) on all admissible d (excluding the degenerate ones)
  rows = select(dd -> !setsearch(Set(DEGEN), dd), ADM);
  KA = matrix(#rows, dv+1, ii, jj,
        if(jj<=dv, VB[jj][rows[ii]+4*r+1], -CD[rows[ii]]));
  my(K2 = matker(KA));
  say(Str("  TEST A (a(d) = kappa c(d) on ", #rows, " admissible d): kernel dim = ", #K2));
  if(#K2==0,
     say("     ==> NO nonzero solution at all: the trace series is NOT proportional to any element of this space."),
     my(haskap = 0);
     for(j=1,#K2, if(K2[dv+1,j]!=0, haskap=j));
     if(haskap==0,
       say(Str("     ==> kernel is entirely kappa=0 (i.e. only forms vanishing on all admissible d): NO match")),
       say(Str("     ==> *** MATCH *** with kappa != 0, column ", haskap, "; kappa = ", K2[dv+1,haskap]))));
  \\ ---- TEST B : square-index dictionary a(m^2) = m beta_s7(m)
  my(NT = min(20, sqrtint(NC)));
  KA = matrix(NT, dv, ii, jj, VB[jj][ii^2+4*r+1]);
  tgt = vectorv(NT, ii, ii*BS7[ii]);
  sol = 0;
  iferr(sol = matinverseimage(KA, tgt), er, sol = 0);
  if(type(sol)=="t_COL" && #sol>0,
     say(Str("  TEST B (a(m^2)=m beta_s7(m), m<=", NT, "): *** EXACT MATCH ***")),
     say(Str("  TEST B (a(m^2)=m beta_s7(m), m<=", NT, "): NO  (rank ", matrank(KA), " of ", NT, "x", dv, ")"));
     kmax=0;
     for(K=1,NT, my(KB=matrix(K,dv,ii,jj,KA[ii,jj]), tb=vectorv(K,ii,tgt[ii]), z=0);
       iferr(z=matinverseimage(KB,tb), e3, z=0);
       if(type(z)=="t_COL" && #z>0, kmax=K, break));
     say(Str("     longest matchable prefix k = ", kmax, "  (k = dim means generic, no hidden solution)")));
  1;
}
run(1, 1);
run(2, 1);
run(1, 0);
run(2, 0);
say("DONE");
quit;
