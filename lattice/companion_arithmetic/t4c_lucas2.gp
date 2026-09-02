default(parisize,"8G");
read("lib.gp");
NT = 900;
IDX = [8, 4];
{ for(u=1,2,
  my(i=IDX[u], R=ROWS[i], r=R[2], AB=genrow(R,NT), A=AB[1], B=AB[2], Nt=vector(NT+1));
  Nt[1]=0; for(n=1,NT, Nt[n+1] = p_dummy);
  print("");
  for(pi=1,4,
    my(p=[5,7,11,13][pi], nmax, ok1=0, ok2=0, ok3=0, tot=0, cand=List(), okc=0, allint=1);
    Nt[1]=0; for(n=1,NT, Nt[n+1] = p^(r*logint(n,p))*B[n+1]);
    for(n=1,NT, if(denominator(Nt[n+1])%p==0, allint=0));
    nmax = floor(NT/p)-1;
    for(n=1,nmax, for(m=0,p-1, tot++;
      if(Mod(Nt[n*p+m+1] - A[n+1]*Nt[m+1] - Nt[n+1]*A[m+1], p)==0, ok1++);
      if(Mod(Nt[n*p+m+1] - A[n+1]*Nt[m+1], p)==0, ok2++);
      if(Mod(Nt[n*p+m+1] - Nt[n+1]*A[m+1], p)==0, ok3++)));
    /* solve for c in  Nt_{np+m} = a_n Nt_m + c Nt_n a_m  using (n,m)=(1,1) */
    my(cc = if(Mod(Nt[2]*A[2],p)==0, -1, lift(Mod(Nt[p+2] - A[2]*Nt[2], p)/Mod(Nt[2]*A[2],p))));
    if(cc>=0, for(n=1,nmax, for(m=0,p-1, if(Mod(Nt[n*p+m+1] - A[n+1]*Nt[m+1] - cc*Nt[n+1]*A[m+1], p)==0, okc++))));
    print("row ", R[1], " p=", p, "  Nt_n := p^(r*floor(log_p n)) b_n  p-integral: ", if(allint,"YES","NO"),
      "  | Nt_{np+m} = a_n Nt_m + Nt_n a_m mod p : ", ok1, "/", tot,
      "  | = a_n Nt_m : ", ok2, "/", tot, "  | = Nt_n a_m : ", ok3, "/", tot,
      "  | = a_n Nt_m + ", cc, "*Nt_n a_m : ", okc, "/", tot)));
}
quit;
