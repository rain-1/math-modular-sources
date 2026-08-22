default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");
\\ Is there ANY 3-adic resource on a Catalan host?  sigma_3 = 0 for every Catalan-class row.
{
my(NN = 400, z, av, bv, x1, x2);
print("=== 3-adic slope on Catalan-class rows (sigma_3 = 0 means no 3-adic resource) ===");
z = row2(12,4,32,NN); av = z[1]; bv = z[2];
print("Zagier E (12,4,32), c = 32 = 2^5:");
foreach([100,200,300,400], n,
  printf("   n=%3d  v2(inc)=%5d (/n = %6.3f)   v3(inc)=%4d\n", n,
     valuation(bv[n+1]/av[n+1] - bv[n]/av[n], 2), valuation(bv[n+1]/av[n+1]-bv[n]/av[n],2)*1.0/n,
     valuation(bv[n+1]/av[n+1] - bv[n]/av[n], 3)));
z = zudrow(300); av = z[1]; bv = z[2];
print("Zudilin Catalan row:");
foreach([100,200,300], n,
  printf("   m=%3d  v2(inc)=%5d (/m = %6.3f)   v3(inc)=%4d\n", n,
     valuation(bv[n+1]/av[n+1] - bv[n]/av[n], 2), valuation(bv[n+1]/av[n+1]-bv[n]/av[n],2)*1.0/n,
     valuation(bv[n+1]/av[n+1] - bv[n]/av[n], 3)));
}
quit;
