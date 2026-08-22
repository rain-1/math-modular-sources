default(parisizemax, 12000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/ap_lib.gp");
NN = 40; MM = 200; mo = 5^MM;
gettime();
r = buildall(NN, MM);
aa = r[1]; bb = r[2]; ee = r[3];
ref = read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/refdata.gp");
ra = ref[1]; rb = ref[2];
okb = 1; oka = 1;
{
for(i = 1, NN+1,
  if(bb[i] != Mod(rb[i], mo), okb = 0; print("b mismatch at n=", i-1));
  if(aa[i] != Mod(ra[i], mo), oka = 0; print("a mismatch at n=", i-1));
);
}
print("checked n=0..", NN);
print("b matches reference: ", okb);
print("a matches reference: ", oka);
print("b_0..b_6 (centered) = ", vector(7, i, centerlift(bb[i])));
print("ref b_0..b_6        = ", vector(7, i, rb[i]));
quit;
