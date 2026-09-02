default(parisize,"8G");
read("lib.gp");
NT = 2000;
PL = [2,3,5,7,11,13,17,19,23,29,31,37,41,43];
{ for(i=13,15,
  my(R=ROWS[i], AB=genrow(R,NT), A=AB[1], B=AB[2], line="");
  for(pi=1,#PL,
    my(p=PL[pi], Nt=vector(NT+1), tot=0, nmax, best=-9, bok=0, cnt, uniq=0);
    Nt[1]=0;
    for(n=1,NT, Nt[n+1] = p^(2*logint(n,p))*B[n+1]);
    nmax = floor(NT/p)-1;
    for(uu=0,p-1,
      cnt=0; tot=0;
      for(n=1,nmax, for(m=0,p-1, tot++; if(Mod(Nt[n*p+m+1] - uu*Nt[n+1]*A[m+1], p)==0, cnt++)));
      if(cnt==tot, uniq++);
      if(cnt>bok, bok=cnt; best=uu));
    line = Str(line, "  ", p, ":", if(bok==tot, Str("u=", lift(Mod(best,p))), Str("FAIL(", bok, "/", tot, ")"))));
  print(R[1], " (k=2)", line));
}
quit;
