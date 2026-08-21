\\ For each spaces/sp_N_w.txt produce prep/pr_N_w.txt:
\\   line 1: N w DIM RANK   (RANK = rank of the a_0=0 sublattice, -1 if no F with a_0=1)
\\   line 2: F0            (integral form with constant term 1, size-reduced)
\\   lines : basis of the a_0=0 sublattice, LLL-reduced w.r.t. coefficients 1..RED, by norm
RED = 40;

readsp(fn) = { my(L, h, dim, M);
  L = readstr(fn);
  h = eval(Str("[", strjoin(strsplit(L[1], " "), ","), "]"));
  dim = h[10];
  if(dim == 0, return([h[2], h[4], h[6], 0, 0]));
  M = matrix(dim, #eval(L[2]), i, j, eval(L[i+1])[j]);
  [h[2], h[4], h[6], dim, M];
}

doprep(fn, out) = { my(R, N, w, ch, dim, M, v, g, HU, U, c, F0, K, KB, T, red, nrm, ord, d);
  R = readsp(fn); N = R[1]; w = R[2]; ch = R[3]; dim = R[4];
  if(dim == 0, write(out, Str(N, " ", w, " ", ch, " 0 -1")); return(0));
  M = R[5];
  v = vectorv(dim, i, M[i,1]);          \\ constant terms
  g = content(v);
  if(g != 1, write(out, Str(N, " ", w, " ", ch, " ", dim, " -1")); return(0));
  HU = mathnf(v~, 1);                    \\ v~ is 1 x dim ; v~ * U = H
  U = HU[2];
  \\ find the column of U giving c with c.v = 1 : HNF puts it last
  c = U[, dim];
  if(c~ * v != 1, c = -c);
  if(c~ * v != 1, write(out, Str(N, " ", w, " ", ch, " ", dim, " -2")); return(0));
  K = matrix(dim, dim-1, i, j, U[i, j]);  \\ remaining columns span the kernel
  d = dim - 1;
  F0 = (c~ * M);
  if(d == 0,
    write(out, Str(N, " ", w, " ", ch, " ", dim, " 0"));
    write(out, Vec(F0));
    return(0));
  KB = K~ * M;                             \\ d x ncoef, rows = kernel forms
  red = matrix(d, RED, i, j, KB[i, j+1]);  \\ truncate to coefficients 1..RED
  T = qflll(red~);                          \\ columns of red~ are the vectors
  KB = (T~) * KB;
  \\ sort rows by truncated norm
  nrm = vector(d, i, sum(j=2, min(RED+1, matsize(KB)[2]), KB[i,j]^2));
  ord = vecsort(nrm, , 1);
  KB = matrix(d, matsize(KB)[2], i, j, KB[ord[i], j]);
  \\ size-reduce F0 against the reduced basis (greedy, truncated norm)
  for(rep = 1, 3,
    for(i = 1, d,
      my(bn = sum(j=2, RED+1, KB[i,j]^2));
      if(bn == 0, next);
      my(ip = sum(j=2, RED+1, F0[j]*KB[i,j]));
      my(q = round(ip/bn));
      if(q != 0, F0 = F0 - q*KB[i,]);
    );
  );
  write(out, Str(N, " ", w, " ", ch, " ", dim, " ", d));
  write(out, Vec(F0));
  for(i = 1, d, write(out, Vec(KB[i,])));
  d;
}
