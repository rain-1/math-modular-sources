/* a4_fold.gp -- the complex "fold" constant xi_+ of a row with complex-conjugate
   characteristic roots, extracted directly from the recurrence:
        xi_n = (B_{n+1} - rho_- B_n)/(A_{n+1} - rho_- A_n)  ->  xi_+  as  n->oo,
   with an asymptotic expansion in 1/n, killed by Neville extrapolation in x=1/n.
   Validation: row eta = (11,5,125), whose fold value THEOREM_B_EXACT.md computes
   modularly as  xi_eta = (1/2)L(chi5,3) + i (Pi/10) L(chi5,2).                  */
default(parisize, 6000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
default(realprecision, 200);
O184=0; for(i=1,#OPS, if(OPS[i][1]=="184", O184=OPS[i]));
rowR3(a,b,c,N)={my(A=vector(N+1),B=vector(N+1)); A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1, A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
              B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3); [A,B];}
N=3000;
pe=rowR3(11,5,125,N);   Ae=pe[1]; Be=pe[2];
pr=aperyPair(O184[4],N); A4=pr[1]; B4=pr[2];

/* Neville extrapolation to x=0 of data (x_j, y_j) */
neville(xs, ys) =
{ my(m=#xs, T=ys);
  for(k=1,m-1, for(i=1,m-k,
    T[i] = (T[i]*(0-xs[i+k]) - T[i+1]*(0-xs[i]))/(xs[i]-xs[i+k])));
  T[1];
}
foldval(A,B,rho, n0, step, K) =
{ my(xs=vector(K), ys=vector(K));
  for(j=1,K, my(n=n0+(j-1)*step);
    xs[j] = 1.0/n;
    ys[j] = (B[n+2]-rho*B[n+1])/(A[n+2]-rho*A[n+1]));
  neville(xs,ys);
}
L53 = lfun(5,3); L52 = lfun(5,2);
print("L(chi5,3) = ", L53);
print("L(chi5,2) = ", L52);
print("\n=== eta = (11,5,125), rho_- = 11-2I ===");
{ forstep(K=10,40,10,
    my(x = foldval(Ae,Be, 11-2*I, 1000, 50, K));
    print("  K=",K,"  xi_+ = ", x);
    print("        Re xi - (1/2)L(chi5,3) = ", real(x)-L53/2);
    print("        Im xi - (Pi/10)L(chi5,2) = ", imag(x)-Pi*L52/10));
}
print("\n=== AESZ 184, rho_- = 44-8I ===");
{ forstep(K=10,40,10,
    my(x = foldval(A4,B4, 44-8*I, 1000, 50, K));
    print("  K=",K,"  xi_+ = ", x);
    print("        Re xi - (1/4)L(chi5,3) = ", real(x)-L53/4);
    print("        Im xi - (Pi/20)L(chi5,2) = ", imag(x)-Pi*L52/20));
}
quit
