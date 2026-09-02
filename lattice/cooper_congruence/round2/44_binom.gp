/* 44_binom.gp -- Task (1): verify closed binomial forms for s10, s7 against AVEC, n<=40. */
read("40_core.gp");

NMAX = 40;

/* s10:  A_n = sum_k binomial(n,k)^4 */
{ f10(n) = my(s=0); for(k=0,n, s += binomial(n,k)^4); s; }

/* s7 :  A_n = sum_k binomial(n,k)^2 binomial(n+k,k) binomial(2k,n) */
{ f7(n) = my(s=0); for(k=0,n, s += binomial(n,k)^2*binomial(n+k,k)*binomial(2*k,n)); s; }

A10 = AVEC(2,NMAX); /* s10 is row 2 */
A7  = AVEC(1,NMAX); /* s7  is row 1 */

print("=== Task (1): closed-form binomial sums vs AVEC, n=0..", NMAX, " ===");

{ ok10 = 1;
  for(n=0,NMAX, my(lhs=A10[n+1], rhs=f10(n));
     if(lhs != rhs, ok10=0; print("s10 MISMATCH at n=",n,": AVEC=",lhs," formula=",rhs)));
  if(ok10, print("s10: MATCH for n=0..",NMAX), print("s10: MISMATCH found (see above)"));
}

{ ok7 = 1;
  for(n=0,NMAX, my(lhs=A7[n+1], rhs=f7(n));
     if(lhs != rhs, ok7=0; print("s7 MISMATCH at n=",n,": AVEC=",lhs," formula=",rhs)));
  if(ok7, print("s7: MATCH for n=0..",NMAX), print("s7: MISMATCH found (see above)"));
}

print("Task (1) done.");
quit;
