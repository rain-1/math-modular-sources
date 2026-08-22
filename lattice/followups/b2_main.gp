/* b2_main.gp -- AESZ 207: archimedean limit to ~330 digits, 2-adic limit to
   ~15000 digits, characteristic roots, tower Frobenius eigenvalues.          */
default(parisize, 12000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
default(realprecision, 420);
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NN=1500;
gettime();
pr=aperyPair(OP[4],NN); A=pr[1]; B=pr[2];
print("built in ", gettime(), " ms;  digits(A_NN)=", #digits(numerator(A[NN+1])));

print("\n=== characteristic roots ===");
{ my(lc = vector(5,i,polcoef(OP[4][i],4,'X)), P=0);
  for(i=1,5, P += lc[i]*x^(4-i));
  print("  char poly = ", P);
  print("  factor    = ", factor(P));
  print("  roots     = ", polroots(P));
  print("  |roots|   = ", vector(4,i,abs(polroots(P)[i])));
  print("  product   = ", polcoef(P,0), " = ", factor(abs(polcoef(P,0))));
}
print("\n=== archimedean limit ===");
{ my(x=B[NN+1]/A[NN+1], d=1.0*x-1.0*B[NN]/A[NN]);
  xi = 1.0*x;
  print("  xi_N   = ", xi);
  print("  last increment = ", d, "  => ~", -log(abs(d))/log(10), " digits");
  write("b2_xi_arch.txt", xi);
}
print("\n=== 2-adic limit ===");
{ my(P2=15000, x=B[NN+1]/A[NN+1]+O(2^P2), y=B[NN]/A[NN]+O(2^P2));
  print("  Cauchy v_2(xi_N - xi_{N-1}) = ", valuation(x-y,2));
  xi2 = x;
  write("b2_xi_2adic.txt", lift(x));
  print("  v_2(xi_2) = ", valuation(x,2));
  print("  slope sigma_2 = ", valuation(x-y,2)*1.0/NN);
}
print("\n=== tower Frobenius eigenvalues at good primes ===");
{ for(j=1,3, my(p=[3,5,7][j]);
    print("  p=",p);
    for(a=1,2,
      my(s=0, prev=0);
      while(a*p^(s+1) <= NN, s++);
      for(t=1,s,
        my(n1=a*p^(t-1), n2=a*p^t);
        if(A[n1+1]==0 || A[n2+1]==0, next);
        my(r = (B[n2+1]/A[n2+1])/(B[n1+1]/A[n1+1]));
        printf("    a=%d  s=%d->%d  ratio=%s  v_p(ratio)=%d  ratio*p^3=%s\n",
               a, t-1, t, r+O(p^6), valuation(r,p), (r*p^3)+O(p^6)))));
}
quit
