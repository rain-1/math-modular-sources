default(parisizemax, 6000000000);
read("sc_rows.gp");
read("sc_diag.gp");
rowdiag("R1 Apery zeta(3)", R1cf, 2, 400, [2,3,5,7,13]);
rowdiag("R2 Zagier C (10,3,9)", R2cf, 2, 400, [2,3,5,7,13]);
rowdiag("R3 AZ eta (11,5,125)", R3cf, 2, 400, [2,3,5,7,13]);
rowdiag("R4 Sym^3 Zagier E", R4cf, 6, 400, [2,3,5,7,13]);
rowdiag("R5 AESZ 207", R5cf, 4, 400, [2,3,5,7,13]);
