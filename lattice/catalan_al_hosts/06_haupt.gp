read("lattice/catalan_al_hosts/lib2.gp");
NQ = 40;
/* calibration: level 16 hauptmodul x = eta_2 eta_16^2 /(eta_1^2 eta_8) */
print("divisors(16)=",divisors(16));
r16 = [-2,1,0,-1,2];
print("orders at cusps 1,2,4,8,16: ", vector(5,i, sum(j=1,5, ligoz(16,divisors(16)[i],divisors(16)[j])*r16[j])));
print("x = ", hauptmod(16,r16,10));
/* level 8: find hauptmodul with zero at cusp 8 (=infty) and pole at each other cusp */
print("\n--- level 8, divisors ", divisors(8));
{for(p=1,3, my(dv=divisors(8), ordv=vector(4)); ordv[4]=1;
   my(wid = 8/(divisors(8)[p]*gcd(divisors(8)[p],8/divisors(8)[p])));
   ordv[p] = -1;  /* try simple pole */
   my(rr = unitexp(8, ordv));
   print("pole at c=",dv[p],"  r=",rr, "  sum r=",sum(i=1,4,rr[i]));
   my(tt = hauptmod(8,rr,12)); print("   t = ", tt));}
quit;
