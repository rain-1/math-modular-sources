default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/census_util.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/multi_prime/out/census_3_introws.log";
NN = 400;
ns = [50, 100, 200, 300, 399];
pl = primes(15);

doit(nm, rr) = { report(fn, nm, rr[1], rr[2], pl, ns); };

logit(fn, Str("### integral / order-2 / order-3 rows, N=", NN, " ###"));

logit(fn, "--- zeta(3) order-3 rows ---");
doit("Domb  row3(10,4,64)", row3(10,4,64,NN));
doit("T     row3(12,4,16)", row3(12,4,16,NN));
doit("Apery row3(17,5,1)",  row3(17,5,1,NN));

logit(fn, "--- Zagier order-2 rows ---");
doit("A row2(7,2,-8)",  row2(7,2,-8,NN));
doit("B row2(9,3,27)",  row2(9,3,27,NN));
doit("C row2(10,3,9)",  row2(10,3,9,NN));
doit("D row2(11,3,-1)", row2(11,3,-1,NN));
doit("E row2(12,4,32)", row2(12,4,32,NN));
doit("F row2(17,6,72)", row2(17,6,72,NN));
doit("s18 cooper18",    cooper18(NN));

logit(fn, "--- AZ order-3 rows ---");
doit("delta row3(7,3,81)",   row3(7,3,81,NN));
doit("zeta  row3(9,3,-27)",  row3(9,3,-27,NN));
doit("eta   row3(11,5,125)", row3(11,5,125,NN));

logit(fn, "DONE");
quit;
