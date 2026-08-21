/* growth.gp -- hypothesis (d) of Theorem F:  v_p(a_n) for the census rows.
   Prints, for each row/slope-prime:
     - the profile v_p(a_n), n<=N
     - max, argmax, and the fit  v_p(a_n) <= C*log_p(n)+C'
     - the minimal value along n, and the subsequence attaining v_p=min
     - comparison with s_p(n) (digit sum) and with v_p(binomial(2n,n))=s_p(n) type laws
   Run:  gp -q -s 4G lattice/theoremF_hyp/growth.gp                              */
default(parisizemax, 8000000000);

if(type(Nlim)!="t_INT", Nlim=3000); N = Nlim;

sdig(n,p) = my(v=digits(n,p)); vecsum(v);

profile(A, p, name, N) =
{ my(v=vector(N,n, valuation(A[n+1],p)), mx=0, am=0, mn=10^9, zeros=0, L);
  for(n=1,N, if(v[n]>mx, mx=v[n]; am=n); if(v[n]<mn, mn=v[n]));
  for(n=1,N, if(v[n]==mn, zeros++));
  print("=== ", name, "  p=",p,"  N=",N);
  print("  min v_p = ", mn, "   attained ", zeros, " times;  first 12 such n: ",
        select(n->v[n]==mn, vector(N,n,n))[1..min(12,zeros)]);
  print("  max v_p = ", mx, " at n=", am, ";   log_p(N) = ", log(N)/log(p));
  /* best constant C in v <= C*log_p(n) */
  my(C=0,nC=0); for(n=2,N, my(r=v[n]*log(p)/log(n)); if(r>C, C=r; nC=n));
  print("  max_n v_p(a_n)/log_p(n) = ", C, " at n=", nC);
  /* compare with digit sum and with number of digits */
  my(dmax=0,ndm=0); for(n=1,N, my(d=v[n]-sdig(n,p)); if(d>dmax,dmax=d;ndm=n));
  print("  max (v_p - s_p(n)) = ", dmax, " at n=", ndm);
  my(emax=0); for(n=1,N, my(d=v[n]-#digits(n,p)); if(d>emax,emax=d));
  print("  max (v_p - #digits_p(n)) = ", emax);
  print("  v_p(a_n), n=1..40: ", v[1..min(40,N)]);
  print("  v_p(a_{p^k}), k=1..: ", vector(floor(log(N)/log(p)),k, v[p^k]));
  print("  v_p(a_{p^k-1}): ", vector(floor(log(N+1)/log(p)),k, if(p^k-1<=N && p^k-1>=1, v[p^k-1], -1)));
  v;
}
