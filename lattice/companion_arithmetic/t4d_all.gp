default(parisize,"8G");
read("lib.gp");
NT = 900;
PL = [2,3,5,7,11,13];
{ for(i=1,#ROWS,
  my(R=ROWS[i], r=R[2], AB=genrow(R,NT), A=AB[1], B=AB[2], line="");
  for(pi=1,#PL,
    my(p=PL[pi], Nt=vector(NT+1), okA=0, okB=0, tot=0, nmax, dmin=99, allint=1);
    Nt[1]=0; for(n=1,NT, Nt[n+1] = p^(r*logint(n,p))*B[n+1]; if(denominator(Nt[n+1])%p==0, allint=0));
    nmax = floor(NT/p)-1;
    for(n=1,nmax, for(m=0,p-1, tot++;
      if(Mod(A[n*p+m+1] - A[n+1]*A[m+1], p)==0, okA++);
      my(dd = Nt[n*p+m+1] - Nt[n+1]*A[m+1]);
      if(Mod(dd,p)==0, okB++);
      if(dd!=0, my(vd=valuation(dd,p)); if(vd<dmin, dmin=vd))));
    line = Str(line, "  p=", p, ": a ", okA, "/", tot, ", b ", okB, "/", tot, " (int:", if(allint,"Y","N"), ",minv:", dmin, ")"));
  print(R[1], " r=", r, " c=", R[5], line));
}
quit;
