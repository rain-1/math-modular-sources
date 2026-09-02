/* 44_s18search.gp -- Task (2): systematic binomial-sum search for s18.
   Building blocks (functions of n,k):
     B1=C(n,k) B2=C(2k,k) B3=C(n+k,k) B4=C(2k,n) B5=C(2n-2k,n-k) B6=C(3k,k) B7=C(n,3k) B8=C(n-k,k)
   extra factor options (index 0..10):
     0: none
     1: 2^k     2: 3^k     3: 4^k
     4: 2^(n-k) 5: 3^(n-k) 6: 4^(n-k)
     7: 2^(n-3k) 8: 3^(n-3k) 9: 4^(n-3k)
     10: (-1)^k
   Exponents e1..e8 in {0,1,2,3}; full brute force over 4^8 * 11 formulas,
   filtered against target A_0..A_6 = 1,6,54,564,6390,76356,948276.        */

NCHK = 6;      /* check against given data n=0..6 */
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

/* opt 0..19: opt 0..9 = plain power factors (or none), opt 10..19 = same times (-1)^k */
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

survivors = List();
count = 0;
tstart = getabstime();

{
EMAX=4;
for(e1=0,EMAX,for(e2=0,EMAX,for(e3=0,EMAX,for(e4=0,EMAX,for(e5=0,EMAX,for(e6=0,EMAX,for(e7=0,EMAX,for(e8=0,EMAX,
  my(e=[e1,e2,e3,e4,e5,e6,e7,e8]);
  count++;
  /* opt-independent product table pw[n+1][k+1] */
  my(pw = matrix(NCHK+1,NCHK+1));
  for(n=0,NCHK, for(k=0,n, my(t=1);
     if(e1>0,t*=Bgrid[1][n+1,k+1]^e1); if(e2>0,t*=Bgrid[2][n+1,k+1]^e2);
     if(e3>0,t*=Bgrid[3][n+1,k+1]^e3); if(e4>0,t*=Bgrid[4][n+1,k+1]^e4);
     if(e5>0,t*=Bgrid[5][n+1,k+1]^e5); if(e6>0,t*=Bgrid[6][n+1,k+1]^e6);
     if(e7>0,t*=Bgrid[7][n+1,k+1]^e7); if(e8>0,t*=Bgrid[8][n+1,k+1]^e8);
     pw[n+1,k+1]=t));
  for(opt=0,19,
    my(ok=1);
    for(n=1,NCHK,
      my(s=0); for(k=0,n, s+=pw[n+1,k+1]*Egrid[opt+1][n+1,k+1]);
      if(s != target[n+1], ok=0; break)
    );
    if(ok, listput(survivors, [e,opt]))
  );
))))))));
}

print("candidates tested (exponent vectors): ", count);
print("survivors matching n=0..",NCHK,": ", #survivors);
for(s=1,#survivors, print("  survivor: e=",survivors[s][1]," opt=",survivors[s][2]));
print("elapsed ms: ", getabstime()-tstart);

quit;
