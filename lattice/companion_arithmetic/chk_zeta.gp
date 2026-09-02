read("lib.gp");
{ my(R=["zetaplus",3,9,3,27,0], AB=genrow(R,40), A=AB[1], bad=0);
  for(n=0,40, if(denominator(A[n+1])!=1, bad=n; break));
  print("(9,3,+27) R3: a_0..a_6 = ", vector(7,j,A[j]), " ; integral to n=40: ", if(bad==0,"YES",Str("NO, first non-integral at n=",bad)));
}
quit;
