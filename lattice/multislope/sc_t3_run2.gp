default(parisizemax, 6000000000);
read("sc_rows.gp");
read("sc_t3.gp");
t3row("R1 Apery zeta(3) [CONTROL]", R1cf, 2, 6000, [2,3,5,7,13], [1,2,3]);
t3row("R4 Sym^3 Zagier E", R4cf, 6, 6000, [2,3,5,7,13], [1,2,3]);
