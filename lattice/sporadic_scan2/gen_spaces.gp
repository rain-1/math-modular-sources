\\ Integral lattice of weight-w forms with RATIONAL q-expansion in the Galois orbit
\\ of one nebentypus character chi mod N (i.e. the Q-rational subspace of
\\ (+)_{sigma} M_w(Gamma_0(N), chi^sigma)), built from Galois traces Tr(zeta^j f)
\\ and saturated to the maximal integral lattice.
\\ Only orbits whose Q-dimension is safely below the q-expansion precision are kept.

charorbits(N, w) = { my(G, L, seen, o, orb, m2);
  G = znstar(N,1); L = List(); seen = vectorsmall(N);
  for(m = 1, N,
    if(gcd(m,N) != 1 || seen[m], next);
    o = znorder(Mod(m,N));
    for(j = 1, o, if(gcd(j,o)==1, m2 = lift(Mod(m,N)^j); seen[m2] = 1));
    if(chareval(G, znconreylog(G,m), -1) == if(w%2==0, 0, 1/2), listput(L, [m, o]));
  );
  Vec(L);
}

orbgens(N, w, m) = { my(G, chi, mf, C, d, K, vy, ph, row, ncoef, gens);
  G = znstar(N,1); chi = znconreylog(G, m);
  mf = mfinit([N, w, [G,chi]], 4);
  if(mfdim(mf) == 0, return([]));
  ncoef = PREC+1;
  C = mfcoefs(mf, PREC);
  d = matsize(C)[2];
  K = mfparams(mf)[5];
  ph = if(type(K)=="t_POL", poldegree(K), 1);
  gens = List();
  if(ph <= 1,
    for(c = 1, d, listput(gens, vector(ncoef, j, simplify(lift(C[j,c])))))
  ,
    vy = variable(K);
    for(c = 1, d,
      for(jj = 0, ph-1,
        row = vector(ncoef, j, my(x = C[j,c]);
                     x = if(type(x)=="t_POLMOD", x, Mod(x, K));
                     simplify(trace(Mod(vy,K)^jj * x)) );
        listput(gens, row)))
  );
  Vec(gens);
}

saturate(gens) = { my(M, H, Y, S);
  if(#gens == 0, return([]));
  M = Mat(vector(#gens, i, gens[i]~));
  M = M * denominator(M);
  M = matrix(matsize(M)[1], matsize(M)[2], i, j, simplify(lift(M[i,j])));
  H = mathnf(M);
  if(matsize(H)[2] == 0, return([]));
  Y = matkerint(H~);
  if(matsize(Y)[2] == 0, S = matid(matsize(H)[1]), S = matkerint(Y~));
  S;
}

\\ saturate using only the first SATN+1 coefficients (enough by the Sturm bound),
\\ then transport the resulting basis to the full PREC+1 coefficients.
saturate_fast(gens, SATN) = { my(ns, S, Mfull, Msmall, dim, out, x);
  if(#gens == 0, return([]));
  ns = min(SATN+1, #gens[1]);
  S = saturate(vector(#gens, i, vector(ns, j, gens[i][j])));
  dim = matsize(S)[2];
  if(dim == 0, return([]));
  Mfull  = Mat(vector(#gens, i, gens[i]~));
  Msmall = matrix(ns, #gens, i, j, Mfull[i,j]);
  out = matrix(#gens[1], dim);
  for(c = 1, dim,
    x = matinverseimage(Msmall, S[,c]);
    if(#x == 0, return(saturate(gens)));
    out[,c] = Mfull * x;
  );
  out;
}
