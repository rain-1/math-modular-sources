/* Independent identification of the p-adic limit eta = -lim a_n/b_n.
   Kubota-Leopoldt:  zeta_p(1+2n) = lim  zeta(1-2K) = lim -B_{2K}/(2K)
   over 2K -> -2n in weight space, i.e. 2K == -2n mod (p-1)p^N, 2K -> oo.
   Claim to test:  zeta_p(2k+1) = 2*eta   with eta = -a_m/b_m (m large).      */
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");
{
kl(p, n, N) =   /* approximation to zeta_p(1+2n) using 2K == -2n mod (p-1)p^N */
  my(md, tk);
  md = (p-1)*p^N;
  tk = md - 2*n;          /* 2K = (p-1)p^N - 2n  >0, == -2n mod md */
  while(tk % 2, tk += md);
  -bernfrac(tk)/tk;
}
{
chk(p, k, m) =
  my(r, et, zz, i);
  r = run0p(p, k, m);
  et = -r[1][m+1]/r[2][m+1];
  print("X_0(", p, ") k=", k, ":  eta = -a_", m, "/b_", m, ";  2*eta should be zeta_", p, "(", 2*k+1, ")");
  for(N = 2, 6,
    zz = kl(p, k, N);
    print("      N=", N, "  (2K=", (p-1)*p^N - 2*k + (((p-1)*p^N-2*k)%2)*(p-1)*p^N, ")   v_", p, "(2*eta - zeta_p) = ",
          if(2*et - zz == 0, "exact", valuation(2*et - zz, p))));
  print("      v_p(2*eta + zeta_p) [wrong-sign control, N=5] = ",
        if(2*et + kl(p,k,5) == 0, "exact", valuation(2*et + kl(p,k,5), p)));
  print("");
}
chk(2,1,60); chk(2,2,60); chk(2,3,60);
chk(3,1,60); chk(3,2,60);
chk(5,1,60); chk(5,2,60);
chk(7,1,60);
quit;
