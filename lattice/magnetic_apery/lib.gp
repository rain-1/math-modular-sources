\\ magnetic_apery: shared library.  (parisizemax must be set by the caller, NOT here:
\\ a default() inside a read() aborts the rest of the file.)
\\ Fricke hosts (hostscan REPORT.md Sec.4.4): [N, C, B, dv, r, tag]
HOSTS = [[5,125,22,[1,5],[-6,6],"N5C125B22"], [6,81,14,[1,2,3,6],[-4,-4,4,4],"N6C81B14"], [6,72,17,[1,2,3,6],[-5,1,-1,5],"N6C72B17_Apery"], [6,64,20,[1,2,3,6],[-6,6,-6,6],"N6C64B20"], [7,49,13,[1,7],[-4,4],"N7C49B13_Cooper_s7"], [8,32,12,[1,2,4,8],[-4,2,-2,4],"N8C32B12_AZeps"], [8,16,24,[1,2,4,8],[-8,16,-16,8],"N8C16B24_Catalan"], [9,27,9,[1,3,9],[-3,0,3],"N9C27B9_AZzeta"], [10,25,6,[1,2,5,10],[-2,-2,2,2],"N10C25B6_Cooper_s10"], [12,9,10,[1,2,3,4,6,12],[-4,4,4,-4,-4,4],"N12C9B10_Domb"], [12,1,34,[1,2,3,4,6,12],[-12,24,12,-12,-24,12],"N12C1B34"], [18,1,14,[1,2,3,6,9,18],[-6,6,12,-12,-6,6],"N18C1B14_Cooper_s18"]];
useries(dv, r, nq) = q*prod(t=1, #dv, eta(q^dv[t] + O(q^nq))^r[t]);
Fseries(dv, r, nq) = 1 - sum(t=1, #dv, r[t]*dv[t]*sum(n=1, (nq-1)\dv[t], sigma(n)*q^(dv[t]*n))) + O(q^nq);
{peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;}
{magfail(S, lim) = my(bad=0); for(m=1, lim, my(c=polcoeff(S,m)); if(denominator(c)!=1 || c%m!=0, bad=m; break)); bad;}
{intfail(S, lim) = my(bad=0); for(m=1, lim, if(denominator(polcoeff(S,m))!=1, bad=m; break)); bad;}
Dinv(S, lim) = sum(m=1, lim, polcoeff(S,m)/m*q^m) + O(q^(lim+1));
Dop(S, lim) = sum(m=0, lim, m*polcoeff(S,m)*q^m) + O(q^(lim+1));
\\ ---- the anti-invariant parametrisation ----
\\ w = u + 1/(C u) is the Hauptmodul of X_0(N)+, v = u - 1/(C u) is anti-invariant.
\\ rho = v * P(w)/Q(w), deg Q = D, deg P <= D-2  <=>  rho = C u (C u^2-1) Pn(u)/Qn(u).
{Qnpoly(qv, C) = my(D); D = #qv - 1; sum(j=0, D, qv[j+1]*(C*U^2+1)^j*(C*U)^(D-j));}
{Pnpoly(k, D, C) = C*U*(C*U^2-1)*(C*U^2+1)^k*(C*U)^(D-2-k);}
{xibasis(C, us, F2, qv, LIM) = my(D, Qs, Qi, r, res, Pp);
  D = #qv - 1; r = D - 1;
  if(r < 1, return(0));
  Qs = subst(Qnpoly(qv,C), U, us);
  Qi = 1/Qs;
  res = vector(r);
  for(k=0, r-1, Pp = subst(Pnpoly(k, D, C), U, us); res[k+1] = Dinv(F2*Pp*Qi, LIM));
  res;}
{phibasis(C, us, F2, qv, LIM) = my(D, Qs, Qi, r, res, Pp);
  D = #qv - 1; r = D - 1;
  Qs = subst(Qnpoly(qv,C), U, us);
  Qi = 1/Qs;
  res = vector(r);
  for(k=0, r-1, Pp = subst(Pnpoly(k, D, C), U, us); res[k+1] = F2*Pp*Qi + O(q^(LIM+1)));
  res;}
{xitomat(Xs, LIM) = my(r); r = #Xs; matrix(r, LIM, j, m, polcoeff(Xs[j], m));}
\\ L = { alpha in Q^r : alpha*T in Z^M },  T[j,m] = coeff of q^m in Xi_j
{latmag(T) = my(r, M, S, T0, Bc, g, de, nn, K, cand);
  r = matsize(T)[1]; M = matsize(T)[2];
  S = List(); T0 = matrix(r,0);
  for(m=1, M, if(#S>=r, break); cand = matconcat([T0, T[,m]]); if(matrank(cand)>matsize(T0)[2], T0=cand; listput(S,m)));
  if(matsize(T0)[2] < r, return(0));
  Bc = T0^(-1);
  for(m=1, M,
    g = Bc*T[,m];
    de = 1; for(i=1,r, de = lcm(de, denominator(g[i])));
    if(de==1, next);
    nn = de*g;
    K = matsolvemod(Mat(nn~), de, 0, 1)[2];
    Bc = K~*Bc;
  );
  Bc;}
\\ ---- general (no Fricke symmetry) parametrisation: rho = u*P(u)/Q(u), deg P <= degQ-2 ----
{xibasisU(us, F2, Qc, LIM) = my(Du, Qs, Qi, r, res);
  Du = #Qc - 1; r = Du - 1;
  if(r < 1, return(0));
  Qs = sum(j=0, Du, Qc[j+1]*us^j);
  Qi = 1/Qs;
  res = vector(r);
  for(k=0, r-1, res[k+1] = Dinv(F2*us^(k+1)*Qi, LIM));
  res;}
{phibasisU(us, F2, Qc, LIM) = my(Du, Qs, Qi, r, res);
  Du = #Qc - 1; r = Du - 1;
  Qs = sum(j=0, Du, Qc[j+1]*us^j);
  Qi = 1/Qs;
  res = vector(r);
  for(k=0, r-1, res[k+1] = F2*us^(k+1)*Qi + O(q^(LIM+1)));
  res;}
\\ widest space on a deg-1-Hauptmodul host: rho = u*P(u)/Q(u), deg P <= degQ-1
\\ (a_0(Phi)=0 at the cusp infinity only; a_0 at the cusp 0 may be nonzero)
{xibasisW(us, F2, Qc, LIM) = my(Du, Qs, Qi, r, res);
  Du = #Qc - 1; r = Du;
  Qs = sum(j=0, Du, Qc[j+1]*us^j);
  Qi = 1/Qs;
  res = vector(r);
  for(k=0, r-1, res[k+1] = Dinv(F2*us^(k+1)*Qi, LIM));
  res;}
\\ ---- CORRECT detection of a magnetic sublattice (its rank may be < r) ----
\\ L(M) = { alpha : alpha*T[,1..M] in Z^M } is a DECREASING family of full-rank lattices;
\\ magnetic elements are exactly those vectors of L(M) whose norm stays BOUNDED as M grows.
{lllrows(L) = my(r, cc, d, Li, Uu);
  r = matsize(L)[1]; cc = matsize(L)[2];
  d = 1; for(i=1,r, for(j=1,cc, d = lcm(d, denominator(L[i,j]))));
  Li = d*L;
  Uu = qflll(Li~);
  ((Li~*Uu)~)/d;}
{rownorm(v) = my(s); s = 0; for(j=1,#v, s = s + v[j]^2); s;}
\\ [min squared norm over an LLL-reduced basis of L(M), that vector, the whole reduced basis]
{magmin(T, M) = my(r, L, Lr, mn, mi);
  r = matsize(T)[1];
  L = latmag(matrix(r, M, j, m, T[j,m]));
  if(type(L)=="t_INT", return([-1, 0, 0]));
  Lr = lllrows(L);
  mn = -1; mi = 0;
  for(i=1, r, my(nv); nv = rownorm(Lr[i,]); if(mn<0 || nv<mn, mn=nv; mi=i));
  [mn, Lr[mi,], Lr];}
{magver(T, al, M) = my(r, ok); r = matsize(T)[1]; ok = M;
  for(m=1, M, if(denominator(sum(j=1,r, al[j]*T[j,m]))!=1, ok = m-1; break)); ok;}
\\ verdict on a space: [1 if a magnetic element survives to MC, witness, minnorm(MA), minnorm(MB)]
{scanT(T, MA, MB, MC) = my(r1, r2, ok);
  r1 = magmin(T, MA); if(r1[1]<0, return([-1, 0, 0, 0]));
  r2 = magmin(T, MB); if(r2[1]<0, return([-1, 0, 0, 0]));
  if(r2[1] != r1[1], return([0, 0, r1[1], r2[1]]));
  ok = magver(T, r2[2], MC);
  [if(ok>=MC, 1, 0), r2[2], r1[1], r2[1]];}
\\ h(W) monic-or-not integral polynomial -> its u-polynomial  u^d h(C u + 1/u)
{touW(hW, C) = my(d); d = poldegree(hW); sum(j=0, d, polcoeff(hW,j)*(C*U^2+1)^j*U^(d-j));}
\\ all LLL-reduced basis vectors of L(M2) that remain magnetic out to M3
{magsub(T, M2, M3) = my(r, L, Lr, out);
  r = matsize(T)[1];
  L = latmag(matrix(r, M2, j, m, T[j,m]));
  if(type(L)=="t_INT", return([]));
  Lr = lllrows(L);
  out = List();
  for(i=1, r, if(magver(T, Lr[i,], M3) >= M3, listput(out, Lr[i,])));
  Vec(out);}
