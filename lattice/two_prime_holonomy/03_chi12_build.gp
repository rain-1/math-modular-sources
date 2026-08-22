default(parisizemax, 8000000000);
default(realprecision, 60);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
{
print("=== conductor-12 well-poised rows for Catalan: construction & identity check ===");
print("G = ", Catalan);
for(bb = 0, 4,
  for(aa = 1, 4*bb,
    my(r, err);
    r = c12row(aa,bb);
    print("a=",aa," b=",bb,"  c1=",r[3],"  c0=",r[4]);
    print("    Q = ", r[1]);
    print("    P = ", r[2]);
    my(dv = c12direct(aa,bb,r[3],r[4],40000), pv = r[1]*Catalan - r[2]);
    print("    sum_{t<=40000} R = ", dv, "   Q*G-P = ", pv, "   diff = ", abs(dv-pv))));
}
quit;
