\\ 09_indep.gp -- Q(x)-linear independence checks on the Hadamard host.
\\ Rank of { x^j f_i : i<=m, j<=Dg } as power series mod a large prime.
default(parisizemax, 12000000000);
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read(concat(outdir, "rows_raw.gp"));
read(concat(outdir, "rows_full.gp"));
read(concat(outdir, "had.gp"));
read(concat(outdir, "cond.gp"));
N = #HZZ; PP = 2^61 - 1;

\\ theta acts coefficientwise: (theta f)_n = n f_n
th(S) = vector(#S, n, n*S[n]);
ONE = vector(N, n, if(n == 1, 0, 0));   \\ the constant 1 has only the n=0 coefficient
\\ we keep an explicit constant term, so work with the coefficient list c_0..c_N.
\\ Build coefficient vectors of length N+1 (index k = coefficient of x^{k-1}).
mk(S, c0) = vector(N+1, k, if(k == 1, c0, S[k-1]));

rk(fs, Dg) =
{ my(L = List(), v, M0, red);
  for(i = 1, #fs,
    for(j = 0, Dg,
      v = vector(N+1, k, if(k-1 >= j, fs[i][k-j], 0));
      listput(L, v)));
  M0 = matrix(#L, N+1, r, c, Mod(numerator(L[r][c]), PP)/Mod(denominator(L[r][c]), PP));
  [#L, matrank(M0)];
}

W1 = mk(WW, 0); TW = mk(th(WW), 0); T2W = mk(th(th(WW)), 0); T3W = mk(th(th(th(WW))), 0);
CD = mk(CND, 0); TC = mk(th(CND), 0);
ONEv = mk(vector(N, n, 0), 1);
HZZv = mk(HZZ, 0);

print("series length N+1 = ", N+1, ", prime = 2^61-1");
{ my(sets = [[ONEv, W1], [ONEv, W1, TW], [ONEv, W1, TW, T2W], [ONEv, W1, TW, T2W, T3W],
             [ONEv, W1, TW, CD, TC], [ONEv, HZZv, W1, TW]],
     nms = ["{1,W}", "{1,W,thW}", "{1,W,thW,th^2W}", "{1,W,thW,th^2W,th^3W}",
            "{1,W,thW,COND,thCOND}", "{1,A_Z(*)A_N,W,thW}"]);
  for(s = 1, #sets,
    print("  ", nms[s], ":");
    for(Dg = 0, 5,
      my(rr = rk(sets[s], Dg));
      printf("     deg<=%d : rank %d / needed %d %s\n", Dg, rr[2], rr[1],
        if(rr[2] == rr[1], "OK", "*** RELATION ***")))); }
\q
