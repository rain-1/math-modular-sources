default(parisize,"10G");
read("lib.gp");
N = 3000;
{ for(i=1,#ROWS,
  my(R=ROWS[i], AB=genrow(R,N), A=AB[1], bad=0, nz=0);
  for(n=0,N, if(denominator(A[n+1])!=1, bad=n; break));
  for(n=0,N, if(A[n+1]==0, nz++));
  print(R[1], ": a_n integral to n=", N, ": ", if(bad==0,"YES",Str("NO at n=",bad)), " ; #{n: a_n=0} = ", nz));
}
quit;
