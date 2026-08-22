/* shared reporting helpers for the census_* scripts.  No parisizemax here. */

/* v_p(denominator(Q_n)) at checkpoints + last-gap slope estimate */
kapline(qv, p, ns) = {
  my(s = "", k = #ns, r);
  for(i = 1, k, s = Str(s, " n=", ns[i], ":", valuation(denominator(qv[ns[i]+1]), p)));
  r = (valuation(denominator(qv[ns[k]+1]), p) - valuation(denominator(qv[ns[k-1]+1]), p)) * 1.0
      / (ns[k] - ns[k-1]);
  [s, r];
};

/* v_p(P_n/Q_n - P_{n-1}/Q_{n-1}) at checkpoints + last-gap slope estimate */
sigline(qv, pv, p, ns) = {
  my(s = "", k = #ns, r, vv = vector(k));
  for(i = 1, k,
    my(n = ns[i], d = pv[n+1]/qv[n+1] - pv[n]/qv[n]);
    vv[i] = if(d == 0, 10^9, valuation(d, p));
    s = Str(s, " n=", n, ":", vv[i]));
  r = (vv[k] - vv[k-1]) * 1.0 / (ns[k] - ns[k-1]);
  [s, r];
};

logit(fn, s) = { print(s); write(fn, s); };

/* full report for one row */
report(fn, name, qv, pv, plist, ns) = {
  my(nn = #qv - 1, allint = 1);
  logit(fn, Str("=== ", name, "   (N=", nn, ", checkpoints ", ns, ") ==="));
  for(n = 0, nn, if(denominator(qv[n+1]) != 1, allint = 0));
  logit(fn, Str("  Q_n integral for all n<=", nn, " ? ", allint));
  for(i = 1, #plist,
    my(p = plist[i], kk = kapline(qv, p, ns), ss = sigline(qv, pv, p, ns));
    logit(fn, Str("  p=", p, "  v_p(den Q_n):", kk[1], "   kappa~", kk[2]));
    logit(fn, Str("  p=", p, "  v_p(incr)   :", ss[1], "   sigma~", ss[2])));
  logit(fn, "");
};
