/* limits.gp -- direct verification of Theorem B (exact form):
      xi = lim B_n/A_n  =  L(Phi, w+1)   for the nine real-fold rows.
   A_n, B_n are generated exactly over Q from the recurrence; no modular
   input is used on the left-hand side.                                   */
read("common.gp");
default(realprecision, 80);

rowsAB(r,a,b,c,NN) = {
  my(A=vector(NN+1), B=vector(NN+1));
  A[1]=1; A[2]=b; B[1]=0; B[2]=1;
  if(r==2,
    for(n=1,NN-1, A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;
                  B[n+2]=((a*n^2+a*n+b)*B[n+1]-c*n^2*B[n])/(n+1)^2),
    for(n=1,NN-1, A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
                  B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3));
  [A,B];
}

NN = 2600;
print("row        NN     xi_NN - L(Phi,w+1)      L(Phi,w+1)");
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], a=R[3], b=R[4], c=R[5]);
  my(P = if(r==2, 1-a*x+c*x^2, 1-2*a*x+c*x^2), rts=polroots(P));
  if(abs(imag(rts[1]))>1e-20, print(nm,"   [complex fold - skipped]"); next);
  my(AB=rowsAB(r,a,b,c,NN), xi=AB[2][NN+1]/AB[1][NN+1]*1.0);
  my(L=Ltarget(nm), d=xi-L);
  print(nm, "\t", NN, "\t", d, "\t", L);
);
}
quit;
