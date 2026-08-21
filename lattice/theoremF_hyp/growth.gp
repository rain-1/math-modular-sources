/* growth.gp -- hypothesis (d) of Theorem F: p-adic valuations v_p(a_n).
   For each census row and slope prime prints:
     min/max of v_p(a_n) over n<=NN, the n attaining the min (this is what
     (d') needs), the best constant in v_p(a_n) <= C*log_p(n), and the
     comparison with the base-p digit sum s_p(n) and digit count.          */
if(type(NN)!="t_INT", NN=3000);

sdig(n,p) = {vecsum(digits(n,p));}

profile(A, p, name, M) =
{ my(v, mx=0, am=0, mn=10^9, cnt=0, S=List(), C=0.0, nC=0, dmax=-10^9, ndm=0, emax=-10^9);
  v = vector(M, n, valuation(A[n+1], p));
  for(n=1, M, if(v[n]>mx, mx=v[n]; am=n); if(v[n]<mn, mn=v[n]));
  for(n=1, M, if(v[n]==mn, cnt++; if(#S<14, listput(S,n))));
  print("=== ", name, "   p=", p, "   NN=", M);
  print("  min v_p = ", mn, "  attained ", cnt, " times; first: ", Vec(S));
  print("  max v_p = ", mx, " at n=", am, ";  log_p(NN)=", log(M)/log(p));
  for(n=2, M, my(r=v[n]*log(p)/log(n)); if(r>C, C=r; nC=n));
  print("  max v_p(a_n)/log_p(n) = ", C, " at n=", nC);
  for(n=1, M, my(d=v[n]-sdig(n,p)); if(d>dmax, dmax=d; ndm=n));
  print("  max (v_p - s_p(n)) = ", dmax, " at n=", ndm);
  for(n=1, M, my(d=v[n]-#digits(n,p)); if(d>emax, emax=d));
  print("  max (v_p - #digits_p(n)) = ", emax);
  print("  v_p(a_n), n=1..40: ", vector(min(40,M), n, v[n]));
  print("  v_p(a_{p^k}), k>=1: ", vector(floor(log(M)/log(p)), k, v[p^k]));
  print("  v_p(a_{p^k-1}), k>=1: ", vector(floor(log(M+1)/log(p)), k, if(p^k-1>=1, v[p^k-1], -1)));
  v;
}
