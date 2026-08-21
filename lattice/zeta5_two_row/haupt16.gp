/* Find eta-quotient Hauptmoduls of Gamma_0(16) (and 8, 32) via Ligozat divisors. */
\p 60
N = 16; D = divisors(N);
/* Ligozat order at cusp a/c (c | N), in the local uniformiser */
ordc(r, c) = (N/(24*gcd(c^2,N))) * sum(i=1,#D, gcd(c,D[i])^2*r[i]/D[i]);
ncusps(N) = sum(i=1,#divisors(N), my(d=divisors(N)[i]); eulerphi(gcd(d,N/d)));
valid(r) = {
  my(s=sum(i=1,#D,r[i]));
  if(s!=0, return(0));                               /* weight 0 */
  if(sum(i=1,#D, D[i]*r[i]) % 24 != 0, return(0));
  if(sum(i=1,#D, (N/D[i])*r[i]) % 24 != 0, return(0));
  my(pr = prod(i=1,#D, D[i]^r[i]));                  /* must be a rational square */
  if(!issquare(pr) && !issquare(1/pr), return(0));
  1; };
print("divisors of ",N,": ",D, "   #cusps=",ncusps(N));
B = 8;
found = List();
{ forvec(r = vector(#D, i, [-B,B]),
  if(valid(r),
    my(o = vector(#D, i, ordc(r, D[i])));
    if(o[#D]==1,
      my(neg = sum(i=1,#D, if(o[i]<0, -o[i], 0)));
      if(neg==1, listput(found, [r,o]))))); }
print("Gamma_0(16) degree-1 eta quotients with a simple zero at oo: ", #found);
for(i=1,min(#found,25), print("  r=",found[i][1],"  cusp orders=",found[i][2]));
quit;
