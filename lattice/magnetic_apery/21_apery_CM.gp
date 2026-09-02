default(parisizemax,12000000000);
read("lib.gp");
NQ = 405; MA=60; MB=150; MC=400;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
print("### Apery host: CM-guided scan.  W-minimal polynomials from 14_Wpoly.out.");
CMN = ["D-8","D-12","D-15","D-20","D-23","D-24fold","D-32","D-36","D-44","D-48","D-72","D-80","cusp"];
{
CMW = [W+16, W+18, W^2+27*W+171, W^2+16*W-16,
 W^6+79*W^5+2229*W^4+20643*W^3-183989*W^2-4820475*W-25493369,
 W^2-288, W^4+68*W^3+1932*W^2+26384*W+140744, W^2+48*W+528,
 W^6+44*W^5+892*W^4+32032*W^3+997296*W^2+13810368*W+67766336,
 W^2-18*W-594, W^2-64*W-1376,
 W^8+216*W^7+27688*W^6+1740320*W^5+59541776*W^4+1187984000*W^3+13920353792*W^2+89323258112*W+243431780416,
 W+17];
}
{
tryh(hW, mult, wide, lab) = my(Qp, qv, Xs, T, res);
  Qp = subst(touW(hW, C), U, x)^mult;
  qv = vector(poldegree(Qp)+1, j, polcoeff(Qp, j-1));
  if(wide, Xs = xibasisW(us, F2, qv, MC), Xs = xibasisU(us, F2, qv, MC));
  if(type(Xs)=="t_INT", print("   ", lab, " dim 0"); return(0));
  T = xitomat(Xs, MC);
  res = scanT(T, MA, MB, MC);
  if(res[1]==1, print("  *** MAGNETIC *** ", lab, "  dim=", #Xs, "  alpha=", res[2]); return(1));
  print("   ", lab, "  dim=", #Xs, "  minnorm ", res[3], " -> ", res[4], "   MISS");
  0;
}
print("-- single CM orbit, pole orders 1,2,3, anti-invariant and widest --");
{
for(i=1, #CMW, for(mm=1, 3,
  tryh(CMW[i], mm, 0, concat(concat("[anti] ", CMN[i]), concat("^", mm)));
  tryh(CMW[i], mm, 1, concat(concat("[wide] ", CMN[i]), concat("^", mm)));
));
}
print("-- pairs, double poles --");
{
for(i=1, #CMW, for(j=i+1, #CMW,
  if(poldegree(CMW[i])+poldegree(CMW[j]) > 4, next);
  tryh(CMW[i]*CMW[j], 2, 0, concat(concat("[anti] ", CMN[i]), concat(" * ", CMN[j])));
  tryh(CMW[i]*CMW[j], 2, 1, concat(concat("[wide] ", CMN[i]), concat(" * ", CMN[j])));
));
}
print("-- triples and quadruples of the small orbits, double poles --");
{
for(i=1, #CMW, for(j=i+1, #CMW, for(k=j+1, #CMW,
  if(poldegree(CMW[i])+poldegree(CMW[j])+poldegree(CMW[k]) > 5, next);
  tryh(CMW[i]*CMW[j]*CMW[k], 2, 0, concat(concat(concat("[anti] ", CMN[i]), concat(" * ", CMN[j])), concat(" * ", CMN[k])));
  tryh(CMW[i]*CMW[j]*CMW[k], 2, 1, concat(concat(concat("[wide] ", CMN[i]), concat(" * ", CMN[j])), concat(" * ", CMN[k])));
)));
}
print("-- the full product of every usable CM orbit of degree <= 2, double poles --");
{
my(P);
P = CMW[1]*CMW[2]*CMW[3]*CMW[4]*CMW[6]*CMW[8]*CMW[10]*CMW[11]*CMW[13];
print("   h(W) = ", P);
tryh(P, 2, 0, "[anti] all deg<=2 CM orbits + cusp, doubled");
tryh(P, 2, 1, "[wide] all deg<=2 CM orbits + cusp, doubled");
}
quit;
