default(parisize,"8G");
read("lib.gp");
NT = 1500;
PL = [2,3,5,7,11,13];
chi1(n) = 1;
chim3(n) = kronecker(-3,n);
chim4(n) = kronecker(-4,n);
chi5(n) = kronecker(5,n);
{ PSI = [chi1, chim3, chim3, chi1, chim4, chim3, chi1, chi1, chi1, chi1, chim3, chi5, chi1, chi1, chim3]; }
{ PSNAME = ["1", "chi_-3", "chi_-3", "1", "chi_-4", "chi_-3", "1", "1", "1", "1", "chi_-3", "chi_5", "?1", "?1", "?chi_-3"]; }
{ for(i=1,#ROWS,
  my(R=ROWS[i], r=R[2], AB=genrow(R,NT), A=AB[1], B=AB[2], ps=PSI[i], line="");
  for(pi=1,#PL,
    my(p=PL[pi], Nt=vector(NT+1), tot=0, nmax, best=-9, bok=0, cnt, tag);
    Nt[1]=0;
    for(n=1,NT, Nt[n+1] = p^(r*logint(n,p))*B[n+1]);
    nmax = floor(NT/p)-1;
    for(uu=0,p-1,
      cnt=0; tot=0;
      for(n=1,nmax, for(m=0,p-1, tot++; if(Mod(Nt[n*p+m+1] - uu*Nt[n+1]*A[m+1], p)==0, cnt++)));
      if(cnt>bok, bok=cnt; best=uu));
    tag = if(Mod(best - ps(p), p)==0, "=psi(p)", Str(" NE psi(p)=", ps(p)));
    line = Str(line, "  p=", p, ": u=", lift(Mod(best,p)), " ", bok, "/", tot, tag));
  print(R[1], " psi=", PSNAME[i], line));
}
quit;
