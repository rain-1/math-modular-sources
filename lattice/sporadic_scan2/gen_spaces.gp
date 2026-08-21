\\ Integral lattice of weight-w forms on Gamma_1(N) with RATIONAL q-expansion.
\\ Built as Q-span of Galois traces Tr_{Q(zeta)/Q}(zeta^j f_i) over all nebentypus
\\ characters chi with chi(-1)=(-1)^w, then saturated to the maximal integral lattice.
\\ Output: one text file per (N,w) with the basis as integer q-expansion vectors.

PREC = 210;          \\ number of q-coefficients (0..PREC)

\\ Galois orbit representatives of Dirichlet characters mod N with given parity
charlist(N, w) = { my(G, L, seen, res, o, orb, m2);
  G = znstar(N,1); L = List(); seen = vectorsmall(N);
  for(m = 1, N,
    if(gcd(m,N) != 1 || seen[m], next);
    o = znorder(Mod(m,N));
    orb = List();
    for(j = 1, o, if(gcd(j,o)==1, m2 = lift(Mod(m,N)^j); listput(orb, m2); seen[m2] = 1));
    \\ parity
    if(chareval(G, znconreylog(G,m), -1) == if(w%2==0, 0, 1/2),
       listput(L, [m, o]));
  );
  Vec(L);
}

ratgens(N, w) = { my(G, CL, gens, mf, C, d, K, vy, ph, row, ncoef);
  G = znstar(N,1);
  CL = charlist(N, w);
  gens = List();
  ncoef = PREC+1;
  for(i = 1, #CL,
    my(m = CL[i][1], chi = znconreylog(G,m));
    mf = mfinit([N, w, [G,chi]], 4);
    if(mfdim(mf) == 0, next);
    C = mfcoefs(mf, PREC);
    d = matsize(C)[2];
    K = mfparams(mf)[5];
    ph = if(type(K)=="t_POL", poldegree(K), 1);
    if(ph <= 1,
      for(c = 1, d, listput(gens, vector(ncoef, j, simplify(lift(C[j,c])))))
    ,
      vy = variable(K);
      for(c = 1, d,
        for(jj = 0, ph-1,
          row = vector(ncoef, j, my(x = C[j,c]);
                       x = if(type(x)=="t_POLMOD", x, Mod(x, K));
                       simplify(trace(Mod(vy,K)^jj * x)) );
          listput(gens, row)
        )
      )
    );
  );
  Vec(gens);
}

\\ saturate: given rational row-generators (as vectors), return integral basis rows
saturate(gens) = { my(M, H, Y, S);
  if(#gens == 0, return([]));
  M = matconcat(vectorv(#gens, i, 0));   \\ placeholder
  M = Mat(vector(#gens, i, gens[i]~));   \\ columns = generators
  M = M * denominator(M);
  M = matrix(matsize(M)[1], matsize(M)[2], i, j, simplify(lift(M[i,j])));
  H = mathnf(M);
  if(matsize(H)[2] == 0, return([]));
  Y = matkerint(H~);
  if(matsize(Y)[2] == 0, S = matid(matsize(H)[1]), S = matkerint(Y~));
  S;    \\ columns = integral basis vectors in Z^ncoef
}

dowork(N, w, fname) = { my(g, S, f, dim);
  g = ratgens(N, w);
  if(#g == 0, write(fname, "DIM 0"); return(0));
  S = saturate(g);
  dim = matsize(S)[2];
  write(fname, Str("N ", N, " w ", w, " DIM ", dim, " PREC ", PREC));
  for(c = 1, dim, write(fname, Vec(S[,c])));
  dim;
}
