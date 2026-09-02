default(parisize,"10G");
read("lib.gp");
NT = 2000;
PL = [2,3,5,7,11,13,17,19,23,29,31,37,41,43];
chi1(n) = 1;
chim3(n) = kronecker(-3,n);
chim4(n) = kronecker(-4,n);
chi5(n) = kronecker(5,n);
{ PSI = [chi1, chim3, chim3, chi1, chim4, chim3, chi1, chi1, chi1, chi1, chim3, chi5, chi1, chi1, chim3]; }
{ KEXP = [2,2,2,2,2,2,3,3,3,3,3,3,2,2,2]; }
{ for(i=1,#ROWS,
  my(R=ROWS[i], k=KEXP[i], AB=genrow(R,NT), A=AB[1], B=AB[2], ps=PSI[i], line="", lmv="");
  for(pi=1,#PL,
    my(p=PL[pi], Nt=vector(NT+1), ok=1, mv=99, good=(R[5]%p!=0 || i>=13));
    if(!good, next);
    Nt[1]=0;
    for(n=1,NT, Nt[n+1] = p^(k*logint(n,p))*B[n+1]);
    /* iterated digit law */
    for(n=1,NT,
      my(dg=digits(n,p), j=#dg-1, rhs=Mod(ps(p)^j*Nt[dg[1]+1],p));
      for(t=2,#dg, rhs = rhs*A[dg[t]+1]);
      if(Mod(Nt[n+1],p)!=rhs, ok=0; break));
    /* sharpness of the one-step law */
    my(nmax=floor(NT/p)-1);
    for(n=1,nmax, for(m=0,p-1,
      my(dd=Nt[n*p+m+1] - ps(p)*Nt[n+1]*A[m+1]);
      if(dd!=0, my(v=valuation(dd,p)); if(v<mv, mv=v))));
    line = Str(line, " ", p, ":", if(ok,"ok","FAIL"));
    lmv = Str(lmv, " ", p, ":", mv));
  print(R[1], " digit law:", line, "  | min v_p(beta_{np+m}-psi(p)beta_n a_m):", lmv));
}
quit;
