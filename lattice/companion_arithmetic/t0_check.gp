default(parisize,"2G");
read("lib.gp");
N = 60;
{ for(i=1,#ROWS,
  my(R=ROWS[i], AB=genrow(R,N), A=AB[1], B=AB[2], bad=0);
  for(n=0,N, if(denominator(A[n+1])!=1, bad=n; break));
  print(R[1], "  r=", R[2], "  c=", R[5], "  a_0..a_5=", vector(6,j,A[j]), "  b_1..b_4=", vector(4,j,B[j+1]), "  A integral: ", if(bad==0,"yes",Str("NO at n=",bad))));
}
quit;
