default(parisizemax, 6000000000);
read("sc_rows.gp");
read("sc_t3b.gp");
t3b("R4 Sym^3 Zagier E [chi test]", R4cf, 6, 15625, [3,5,7,11,13,17,19,23], [1,2,3]);
