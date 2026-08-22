default(parisizemax, 6000000000);
read("sc_rows.gp");
{ gettime(); a = genseq(R4cf, 6, [1], 8192); print("R4 A N=8192: ", gettime(), " ms, digits ", #digits(a[8193]));
  x = compan(R4cf, 6, 3, 8192); print("R4 X3 N=8192: ", gettime(), " ms, den digits ", #digits(denominator(x[8193]))); }
