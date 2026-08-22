default(parisize, 8000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NN=500; pr=aperyPair(OP[4],NN); A=pr[1]; B=pr[2];
x = B[NN+1]/A[NN+1];
print("v_2 = ", valuation(x,2));
print("16*xi_2 mod 2^120 = ", lift((16*x)+O(2^120)));
print("16*xi_2 mod 2^40  = ", lift((16*x)+O(2^40)));
print("Cauchy at N=500: ", valuation(x - B[NN]/A[NN], 2));
quit
