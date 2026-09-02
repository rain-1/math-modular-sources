default(parisize,"10G");
read("lib.gp");
NT = 2000;
PL = [2,3,5,7,11,13,17,19,23,29,31,37,41,43];
{ KEXP = [2,2,2,2,2,2,3,3,3,3,3,3,2,2,2]; }
iszp(x,p) = (denominator(x)%p != 0);
{ for(i=1,#ROWS,
  my(R=ROWS[i], k=KEXP[i], AB=genrow(R,NT), A=AB[1], B=AB[2], lg="", lb="");
  for(pi=1,#PL,
    my(p=PL[pi], Nt=vector(NT+1), tot=0, nmax, good=(R[5]%p!=0), sols=List(), bok=0);
    Nt[1]=0;
    for(n=1,NT, Nt[n+1] = if(good, p^(k*logint(n,p))*B[n+1], B[n+1]));
    nmax = floor(NT/p)-1;
    for(uu=0,p-1,
      my(cnt=0); tot=0;
      for(n=1,nmax, for(m=0,p-1,
        my(dd = Nt[n*p+m+1] - uu*Nt[n+1]*A[m+1]);
        if(iszp(dd,p), tot++; if(dd==0 || valuation(dd,p)>=1, cnt++))));
      if(cnt==tot, listput(sols,uu));
      if(cnt>bok, bok=cnt));
    my(s = Str(" ", p, ":", if(#sols>0, Str("u=", Vec(sols)), Str("X(", bok, "/", tot, ")"))));
    if(good, lg=Str(lg,s), lb=Str(lb,s)));
  print(R[1], " k=", k, " c=", R[5], " | GOOD p:", lg, " | BAD p (b unnormalised):", lb));
}
quit;
