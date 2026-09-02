\\ 10_flat.gp -- is Phi_flat = sum_e mu(e)psi(e) e Phi(e tau) (coefficients n beta(n) = a(n^2))
\\ a meromorphic modular form on Gamma_0(N)?  Test: is Phi_flat/Phi a rational function of u?
read("lib.gp");
M = 160;
{
for(k=1,3,
  my(S=Setup(k,M+4), u=S[1], Ph=S[4], c, flat, PF, r, mat, ker, DD);
  c = vector(M,m,polcoeff(Ph,m));
  flat = vector(M);
  for(n=1,M, my(s=0); fordiv(n,e, s += moebius(e)*psin(k,e)*e*c[n/e]); flat[n]=s);
  PF = sum(n=1,M, flat[n]*'q^n) + O('q^(M+1));
  r = PF/(Ph+O('q^(M+1)));       \\ q-series of the ratio, starts at 1
  \\ test: r = P(u)/Q(u) with deg <= DD ; i.e. sum_{i<=DD} a_i u^i - r * sum_{j<=DD} b_j u^j = 0
  DD = 6;
  mat = matrix(M-2, 2*DD+2, i, j,
     if(j<=DD+1, polcoeff(u^(j-1)+O('q^(M+1)), i), -polcoeff(r*u^(j-DD-2)+O('q^(M+1)), i)));
  ker = matker(mat);
  print("row ",NAM[k],"  Phi_flat/Phi rational in u of degree <= ",DD,"?  dim ker = ",#ker);
  print("   ratio q-expansion: ", vector(8,i,polcoeff(r,i-1)));
);
}
quit;
