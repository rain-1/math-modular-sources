/* exactlaws.gp -- the two exact digit laws, pushed to n <= 10000.
     B at p=3 :  v_3(a_n) = s_3(n)
     E at p=2 :  v_2(a_n) = 2 s_2(n)                                       */
read("lattice/euler_criterion/rows.gp");
M = 10000;
{ my(A, ok);
  A = rowR2(9,3,27,M)[1]; ok = 1;
  for(n=1,M, if(valuation(A[n+1],3) != vecsum(digits(n,3)), ok=0; print("B FAIL n=",n); break));
  print("B  v_3(a_n) = s_3(n)  for 1<=n<=", M, " : ", if(ok,"TRUE","FALSE"));
  A = rowR2(12,4,32,M)[1]; ok = 1;
  for(n=1,M, if(valuation(A[n+1],2) != 2*vecsum(digits(n,2)), ok=0; print("E FAIL n=",n); break));
  print("E  v_2(a_n) = 2 s_2(n) for 1<=n<=", M, " : ", if(ok,"TRUE","FALSE"));
}
quit;
