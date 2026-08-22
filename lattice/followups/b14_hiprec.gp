default(parisize, 20000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
default(realprecision, 1250);
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NN=5000;
gettime(); pr=aperyPair(OP[4],NN); A=pr[1]; B=pr[2];
print("built ", gettime()," ms; digits A_NN=",#digits(numerator(A[NN+1])));
xi = 1.0*B[NN+1]/A[NN+1];
d  = xi - 1.0*B[NN]/A[NN];
print("increment ~ 10^", log(abs(d))/log(10));
write("b14_xi_arch_1200.txt", xi);
print("xi = ", xi);
quit
