default(parisizemax, 8000000000);
\\ 02_run_all.gp -- consolidated run: exact joint recurrences + singularity data
\\ for the conductor-6 (chi6row) and conductor-12 (c12rowB, and its two halves)
\\ well-poised families at alpha = 2 and alpha = 1.
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/02_fit_recur.gp");
\\ sel: 0 = joint (Q,P), 1 = Q only, 2 = P only
run("c6_a2",     2, 20, 0);
run("c6_a2",     1, 12, 1);
run("c6_a1",     2, 20, 0);
run("c12_a2",    1, 22, 1);
run("h12_a2_r1", 2, 20, 0);
run("h12_a2_r5", 2, 20, 0);
run("h12_a1_r1", 2, 20, 0);
run("h12_a1_r5", 2, 20, 0);
quit;
