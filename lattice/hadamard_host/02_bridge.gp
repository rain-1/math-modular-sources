\\ 02_bridge.gp -- reproduce v_2(a^Z b^N - a^N b^Z) >= 24n - O(log n)
\\ in BOTH normalisations: the published integer rows (X,Y;V,U) and the
\\ rational rows (az,bz;an,bn) that define the Hadamard host.
default(parisizemax, 8000000000);
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read(concat(outdir, "csvrows.gp"));
read(concat(outdir, "rows_raw.gp"));
NC = #BN0;
DT = vector(6*NC+1); DT[1] = 1;
for(k = 1, 6*NC, DT[k+1] = lcm(DT[k], k));
Dn(n) = DT[6*n+1];
ez(n) = my(m = 3*n); min(6*m, 4*m + 3 + logint(2*m-1, 2));

print("n  v2(h)/n   v2(h)-24n   v2(w)/n   v2(w)+2n   check");
HH = vector(NC); WW = vector(NC);
{ for(n = 1, NC,
    my(h = XCSV[n]*UCSV[n] - VCSV[n]*YCSV[n],
       w = AZ[n]*BN0[n] - AN[n]*BZ[n], ok);
    HH[n] = h; WW[n] = w;
    ok = (h == 2^(ez(n)+14*n) * Dn(n)^4 * w);
    if(n % 7 == 0 || n <= 4,
      printf("%3d  %8.4f  %+9.3f  %8.4f  %+9.3f  %s\n",
        n, valuation(h,2)/n, valuation(h,2)-24*n,
        valuation(numerator(w),2)-valuation(denominator(w),2),
        0, if(ok,"ok","BAD")));
  ); }
print("");
print("v_2 profile of the rational Hadamard minor w_n = az_n*bn_n - an_n*bz_n:");
{ for(n = 1, NC, if(n % 7 == 0 || n <= 4,
    printf("  n=%3d  v2(w_n)=%7d   v2(w_n)/n=%9.5f   v2(h_n)=%7d  v2(h_n)/n=%9.5f  ez+14n=%6d  4v2(D)=%3d\n",
      n, valuation(WW[n],2), valuation(WW[n],2)/n,
      valuation(HH[n],2), valuation(HH[n],2)/n, ez(n)+14*n, 4*valuation(Dn(n),2)))); }
print("");
{ my(vs = vector(NC, n, valuation(HH[n],2)), best = 10^9, bn = 0);
  for(n = 1, NC, my(dev = vs[n] - 24*n); if(dev < best, best = dev; bn = n));
  print("min over n<=", NC, " of  v2(h_n) - 24n  = ", best, "  at n = ", bn);
  print("v2(h_n) - 24n at n = 20,40,60,80,98: ",
        vector(5, k, my(n=[20,40,60,80,98][k]); vs[n]-24*n));
  print("least-squares slope of v2(h_n) on n in [50,98]: ");
  my(sx=0,sy=0,sxx=0,sxy=0,cnt=0);
  for(n=50,NC, sx+=n; sy+=vs[n]; sxx+=n^2; sxy+=n*vs[n]; cnt++);
  printf("   slope = %.6f   intercept = %.4f\n",
    (cnt*sxy-sx*sy)/(cnt*sxx-sx^2), (sy*sxx-sx*sxy)/(cnt*sxx-sx^2));
}
write(concat(outdir,"minor.gp"), "WW = ", WW, ";");
\q
