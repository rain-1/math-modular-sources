default(parisizemax, 6000000000);
read("sc_rows.gp");
read("sc_t2c.gp");
t2c("R5 AESZ 207", R5cf, 4, 2, 12, 4000, 40, 2);
t2c("R5 AESZ 207 (s>=6)", R5cf, 4, 2, 12, 4000, 40, 6);
t2c("R2 Zagier C", R2cf, 2, 3, 2, 4000, 40, 2);
t2c("R3 AZ eta", R3cf, 2, 5, 3, 4000, 40, 2);
