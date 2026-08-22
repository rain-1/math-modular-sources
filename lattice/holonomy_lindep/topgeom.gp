default(parisizemax, 6G);
read("/home/ubuntu/code/math-modular-sources/lattice/mum_survey/ops.gp");
wantv = ["243","388","101","313","34","104","380","102","113","103","206","207","410","392"];
wants = Set(wantv);
for(i=1,#OPS, my(o=OPS[i], nm=o[1]); if(setsearch(wants, nm), print("AESZ ", nm, " | cluster ", o[2], " | degz ", o[3], " | sing z: ", o[8], " | A_n: ", vector(min(6,#o[9]),j,o[9][j]))));
