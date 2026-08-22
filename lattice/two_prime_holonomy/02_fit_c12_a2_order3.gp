default(parisizemax, 8000000000);
\\ 02_fit_c12_a2_order3.gp -- the MINIMAL joint (Q,P) recurrence for the
\\ conductor-12 product-form family c12rowB(2b,2b,b): order 3, degree 60,
\\ nullity 1, chi(x) = (x - 1048576)(x + 1024)^2.
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/02_fit_recur.gp");
run("c12_a2", 3, 62, 0);
quit;
