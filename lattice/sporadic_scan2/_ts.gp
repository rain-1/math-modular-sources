default(parisizemax, 4000000000);
read("gen_spaces.gp");
PREC = 600;
t0=getabstime();
g = orbgens(12, 1, 5);
print("#g=",#g, " t=",getabstime()-t0);
t0=getabstime();
S = saturate_fast(g, 210);
print("dim=",matsize(S)[2], " t=",getabstime()-t0, " len=",matsize(S)[1]);
print(Vec(S[,1])[1..10]);
quit
