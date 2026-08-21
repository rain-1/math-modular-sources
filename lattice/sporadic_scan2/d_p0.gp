default(parisizemax, 5000000000);
read("gen_spaces.gp");
read("prep.gp");
PREC = 600;
SPDIR="spaces600";
PRDIR="prep600";
LEVELS = [4, 10, 18];
read("run_prep.gp");
