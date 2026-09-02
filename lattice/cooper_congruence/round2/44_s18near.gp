/* 44_s18near.gp -- Task (2) continued: report best near-misses (longest matching
   prefix n=0..k, or smallest relative miss) over the same search space as
   44_s18search.gp (EMAX=3, 20 extra-factor options), since the exhaustive
   search (44_s18search.gp / EMAX up to 5) found NO exact match on n=0..6.   */

NCHK = 6;
target = [1,6,54,564,6390,76356,948276];

{ Bval(i,n,k) =
    if(i==1, binomial(n,k),
    if(i==2, binomial(2*k,k),
    if(i==3, binomial(n+k,k),
    if(i==4, binomial(2*k,n),
    if(i==5, binomial(2*n-2*k,n-k),
    if(i==6, binomial(3*k,k),
    if(i==7, binomial(n,3*k),
    binomial(n-k,k)))))))); /* i==8 */
}
{ Bgrid = vector(8, i, matrix(NCHK+1,NCHK+1, np,kp, if(kp-1>np-1,0,Bval(i,np-1,kp-1)))); }

{ Extra(opt,n,k) = my(m=opt%10, p);
  p = if(m==0, 1,
      if(m==1, 2^k,
      if(m==2, 3^k,
      if(m==3, 4^k,
      if(m==4, 2^(n-k),
      if(m==5, 3^(n-k),
      if(m==6, 4^(n-k),
      if(m==7, 2^(n-3*k),
      if(m==8, 3^(n-3*k),
      4^(n-3*k) )))))))));
  if(opt>=10, p*(-1)^k, p);
}
{ Egrid = vector(20, o, matrix(NCHK+1,NCHK+1, np,kp, if(kp-1>np-1,0,Extra(o-1,np-1,kp-1)))); }

bestlen = -1;
bestlist = List();

{
EMAX=3;
for(e1=0,EMAX,for(e2=0,EMAX,for(e3=0,EMAX,for(e4=0,EMAX,for(e5=0,EMAX,for(e6=0,EMAX,for(e7=0,EMAX,for(e8=0,EMAX,
  my(e=[e1,e2,e3,e4,e5,e6,e7,e8]);
  my(pw = matrix(NCHK+1,NCHK+1));
  for(n=0,NCHK, for(k=0,n, my(t=1);
     if(e1>0,t*=Bgrid[1][n+1,k+1]^e1); if(e2>0,t*=Bgrid[2][n+1,k+1]^e2);
     if(e3>0,t*=Bgrid[3][n+1,k+1]^e3); if(e4>0,t*=Bgrid[4][n+1,k+1]^e4);
     if(e5>0,t*=Bgrid[5][n+1,k+1]^e5); if(e6>0,t*=Bgrid[6][n+1,k+1]^e6);
     if(e7>0,t*=Bgrid[7][n+1,k+1]^e7); if(e8>0,t*=Bgrid[8][n+1,k+1]^e8);
     pw[n+1,k+1]=t));
  for(opt=0,19,
    my(len=0);
    for(n=1,NCHK,
      my(s=0); for(k=0,n, s+=pw[n+1,k+1]*Egrid[opt+1][n+1,k+1]);
      if(s == target[n+1], len++, break)
    );
    if(len>bestlen, bestlen=len; bestlist=List([[e,opt,len]]),
       if(len==bestlen && len>0, listput(bestlist,[e,opt,len])))
  );
))))))));
}

print("best matching-prefix length (n=1..): ", bestlen, " out of ", NCHK);
print("number of candidates achieving it: ", #bestlist);
kshow = min(#bestlist, 15);
{ for(i=1,kshow,
  my(e=bestlist[i][1], opt=bestlist[i][2]);
  print("  e=",e," opt=",opt);
); }

/* for the top few, print actual A_n produced vs target for n=0..NCHK, plus n=7,8 to see divergence */
NEXT = 8;
{ Bgrid2 = vector(8, i, matrix(NEXT+1,NEXT+1, np,kp, if(kp-1>np-1,0,Bval(i,np-1,kp-1)))); }
{ Egrid2 = vector(20, o, matrix(NEXT+1,NEXT+1, np,kp, if(kp-1>np-1,0,Extra(o-1,np-1,kp-1)))); }
{ evalcand(e,opt,NMAXX) = vector(NMAXX+1,np, my(n=np-1,s=0);
    for(k=0,n, my(t=Egrid2[opt+1][n+1,k+1]);
      for(i=1,8, if(e[i]>0, t*=Bgrid2[i][n+1,k+1]^e[i]));
      s+=t);
    s);
}
print("\n--- sample candidate values vs target (n=0..",NEXT,") ---");
{ for(i=1,kshow,
  my(e=bestlist[i][1], opt=bestlist[i][2]);
  print("e=",e," opt=",opt," -> ", evalcand(e,opt,NEXT));
); }
print("target(n=0..",NCHK,")=",target);

quit;
