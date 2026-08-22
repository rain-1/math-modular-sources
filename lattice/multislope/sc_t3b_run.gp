default(parisizemax, 6000000000);
read("sc_rows.gp");
read("sc_t3b.gp");
t3b("R1 Apery zeta(3) [CONTROL]", R1cf, 2, 8192, [2,3,5,7,13], [1,2,3]);
t3b("R4 Sym^3 Zagier E", R4cf, 6, 8192, [2,3,5,7,13], [1,2,3]);
t3b("R5 AESZ 207", R5cf, 4, 8192, [3,5,7], [1,2,3]);
