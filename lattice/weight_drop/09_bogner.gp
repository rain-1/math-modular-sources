/* (n+1) | A_n for the three Cooper rows, and A_n integral, n<=600 */
{par(ab,N)=my(a=ab[1],b=ab[2],c=ab[3],d=ab[4],A=vector(N+1));A[1]=1;A[2]=b;
 for(n=1,N-1,A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]+n*(c*n^2+d)*A[n])/(n+1)^3);A;}
rows=[["s7",[13,4,27,-3]],["s10",[6,2,64,-4]],["s18",[14,6,-192,12]]];
N=600;
{for(r=1,3, my(A=par(rows[r][2],N));
  print(rows[r][1], ": A_n in Z? ", vecmax(vector(N+1,k,denominator(A[k])))==1,
        "   (n+1)|A_n? ", vecmax(vector(N+1,k,denominator(A[k]/k)))==1,
        "   min v_2(A_n), n>=1 = ", vecmin(vector(N,k,valuation(A[k+1],2))));
);}
quit;
