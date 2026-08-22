/* ap_xcheck.gp -- independent EXACT-RATIONAL recomputation of a_n,b_n for n<=NX
   using the completely different route of cdt_ab.gp (serreverse of x(q) and
   subst), then comparison with the Z/5^MM peeling result.  */
default(parisizemax, 12000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/ap_lib.gp");
NX = 600; MM = 2000; mo = 5^MM;
gettime();
ex = run0p(5, 1, NX, 1);      /* exact rationals, nrm=1 (undivided E_2^*) */
print("exact rational run to n=", NX, "  ms=", gettime());
ap = buildall(NX, MM);
oka = 1; okb = 1; oke = 1;
{
for(i = 1, NX+1,
  if(ap[1][i] != Mod(ex[1][i], mo), oka = 0; print("a mismatch n=", i-1));
  if(ap[2][i] != Mod(ex[2][i], mo), okb = 0; print("b mismatch n=", i-1));
);
}
print("EXACT vs mod-5^", MM, " agreement for 0<=n<=", NX, ":  a:", oka, "  b:", okb);
/* also: exact 5-adic valuations from the exact rationals */
print("v5(b_n) exact, n=590..600: ", vector(11, i, my(z=ex[2][590+i]); if(z==0, -1, valuation(z,5))));
print("v5(b_n) from mod run    : ", vector(11, i, v5(ap[2][590+i], MM)));
quit;
