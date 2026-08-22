read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
pr=aperyPair(OP[4],14); A=pr[1];
for(n=0,10, print("A_",n," = ",A[n+1],"   = ",if(A[n+1]==0,"0",factor(abs(A[n+1])))));
print("\nA_n / (-1)^n:");
for(n=0,8, print("  ", A[n+1]*(-1)^n));
print("\nratios A_n/A_{n-1}:");
for(n=1,10, print("  ", 1.0*A[n+1]/A[n]));
print("\nA_n / binom(2n,n):");
for(n=0,8, print("  ", A[n+1]/binomial(2*n,n), "  integral? ", denominator(A[n+1]/binomial(2*n,n))==1));
print("\nA_n / binom(2n,n)^2:");
for(n=0,8, print("  ", A[n+1]/binomial(2*n,n)^2, "  integral? ", denominator(A[n+1]/binomial(2*n,n)^2)==1));
print("\nleading char poly factorisation of the operator: ");
print(factor(x^4 - 17152*x^3 - 6696206336*x^2 + 255108172480512*x - 47569271064100864));
print("53248=",factor(53248)," ; 89344=",factor(89344)," ; 2^24=16777216");
print("z_-*z_+ = -1/2^24 ; z_-+z_+ = 89344/2^24 = 349/65536");
quit
