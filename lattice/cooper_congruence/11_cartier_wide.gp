/* 11_cartier_wide.gp -- wide-range verification of the mechanism
      C(eta) = psi(p) eta  mod p,   i.e.  a_{pj+p-1} = psi(p) a_j  mod p,
   and of its j=0 slot  a_{p-1} = psi(p) mod p.                                */
default(parisize, 8000000000);
read("lib.gp");
NMAX = 1200;
gettime();
{ AX = vector(3);
  for(k=1,3, my(R=ROWS[k], B=R[3], C=R[4], A=genrow(k,NMAX)[1], F, l, P);
    F = Ser(A,'x);
    l = 'x*Ser(vector(NMAX+1,n,A[n]/n),'x);
    P = 1 - 2*B*'x + (B^2-4*C)*'x^2 + O('x^(NMAX+1));
    AX[k] = l/('x*sqrt(P)*F)); }
print("eta built to x^",NMAX,",  ",gettime()," ms");
av(k,j) = polcoeff(AX[k],j);
print();
print("=== a_{p-1} = psi(p) mod p,  p <= 199 ===");
{ for(k=1,3, my(bad=[]);
   forprime(p=2,199, if((av(k,p-1)-psival(k,p))%p!=0, bad=concat(bad,[p])));
   print("  ",NAM[k],": ",if(#bad==0,"HOLDS for all p<=199",Str("FAILS at ",bad)))); }
print();
print("=== C(eta) = psi(p) eta mod p:  [p, v_p(gcd of defects), #tests],  p<=101 ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,101, my(ps=psival(k,p), g=0, cnt=0);
     for(j=0,(NMAX-p+1)\p, my(i=p*j+p-1); if(i>NMAX,break); cnt++;
       g=gcd(g, av(k,i)-ps*av(k,j)));
     row=concat(row,[[p, if(g==0,"inf",valuation(g,p)), cnt]]));
   print("  ",NAM[k],":");
   print("    ",row)); }
quit;
