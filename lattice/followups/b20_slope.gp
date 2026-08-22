default(parisize, 8000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NN=1600; pr=aperyPair(OP[4],NN); A=pr[1]; B=pr[2];
PR=25000;
print("2-adic Cauchy precision v_2(xi_N - xi_{N-1}) and slope:");
prev=0;
{ forstep(N=200,1600,200,
    my(c=valuation(B[N+1]/A[N+1]-B[N]/A[N],2));
    printf("  N=%d  v_2 = %d  v_2/N = %.4f", N, c, c*1.0/N);
    if(prev, printf("   incr/200 = %.4f", (c-prev)/200.));
    print(); prev=c); }
print("\nv_2(A_n) sample (kappa_2):");
forstep(n=200,1600,200, print("  n=",n,"  v_2(A_n)=",valuation(A[n+1],2)));
print("\nv_13, v_17 slopes (bad primes):");
{ foreach([13,17], p,
    printf("  p=%d: ",p);
    forstep(N=400,1600,400, printf(" v(%d)=%d",N,valuation(B[N+1]/A[N+1]-B[N]/A[N],p)));
    print()); }
quit
