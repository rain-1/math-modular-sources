/* ap_run.gp -- main Z/5^MM run: build a_n,b_n,e_n and dump 5-adic valuations. */
default(parisizemax, 12000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/ap_lib.gp");
NN = 500;
MM = 3*NN + 150;
mo = 5^MM;
outf = concat(concat("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/ap_out_", NN), ".txt");
gettime();
r = buildall(NN, MM);
aa = r[1]; bb = r[2]; ee = r[3];
print("build done, ms=", gettime());

/* --- eta from the top index --- */
vbN = v5(bb[NN+1], MM);
etaN = -Mod(lift(aa[NN+1])/5^vbN, mo) * Mod(lift(bb[NN+1])/5^vbN, mo)^(-1);
print("v5(b_N) = ", vbN, "   (eta determined mod 5^", MM-vbN, ")");
print("eta_N mod 5^40 = ", lift(etaN) % 5^40);
print("eta_N 5-adic digits (first 30): ", digits(lift(etaN) % 5^30, 5));

write(outf, "# NN=", NN, " MM=", MM);
write(outf, "# eta_N mod 5^", MM-vbN, " = ", lift(etaN) % 5^(MM-vbN));
write(outf, "# n  v5(a_n)  v5(b_n)  v5(e_n)  v5(h_n)");
{
for(n = 0, NN,
  my(hn, va, vb, ve, vh);
  hn = aa[n+1] + etaN*bb[n+1];
  va = v5(aa[n+1], MM); vb = v5(bb[n+1], MM); ve = v5(ee[n+1], MM); vh = v5(hn, MM);
  write(outf, n, " ", va, " ", vb, " ", ve, " ", vh);
);
}
print("dump done, ms=", gettime());
/* also save the raw top coefficients for the KL cross-check */
write(outf, "# a_N = ", lift(aa[NN+1]));
write(outf, "# b_N = ", lift(bb[NN+1]));
quit;
