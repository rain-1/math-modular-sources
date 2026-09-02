default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 60);
NQ = 262; NA = 250;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
xs = us/(1+B*us+C*us^2);
Av = peel2(Fs, xs, NA, NQ);
print("### THE UNIT OBSTRUCTION");
print("");
print("(1) identity  x(tau_0) = -1/(W_0 + B)  for a pole of  1 + B' u + C u^2,  W_0 = C u + 1/u = -B'");
print("    check symbolically:  1 + B u + C u^2  mod (C u^2 + B' u + 1)  =  ", lift(Mod(1+B*x+C*x^2, C*x^2+Bp*x+1)));
print("    so 1+Bu+Cu^2 = (B-B')u on the polar locus and x = u/((B-B')u) = 1/(B-B') = -1/(W_0+B).");
print("");
print("(2) |x(tau_0)| over the Galois orbit of each CM point of X_0(6), B = 17:");
print("    minpoly of Y = W+17 ; norm N(Y) = +-1 or an integer; |x_i| = 1/|Y_i|");
CMN = ["D-8","D-12","D-15","D-20","D-23","D-24 fold","D-32","D-36","D-44","D-48","D-72","D-80","D-96","D-120"];
{
CMW = [W+16, W+18, W^2+27*W+171, W^2+16*W-16,
 W^6+79*W^5+2229*W^4+20643*W^3-183989*W^2-4820475*W-25493369,
 W^2-288, W^4+68*W^3+1932*W^2+26384*W+140744, W^2+48*W+528,
 W^6+44*W^5+892*W^4+32032*W^3+997296*W^2+13810368*W+67766336,
 W^2-18*W-594, W^2-64*W-1376,
 W^8+216*W^7+27688*W^6+1740320*W^5+59541776*W^4+1187984000*W^3+13920353792*W^2+89323258112*W+243431780416,
 W^4-108*W^3-8244*W^2-167184*W-1073736, W^2-288*W-5184];
}
print("");
print("disc | minpoly of Y=W+17 | N(Y) | |Y_i| (=1/|x_i|) | max |Y_i| | tail rate");
{
for(i=1, #CMW,
  my(PY, rts, mx, nrm);
  PY = subst(CMW[i], W, Y-17);
  rts = polroots(PY);
  nrm = polcoeff(PY,0)*(-1)^poldegree(PY);
  mx = 0; for(j=1,#rts, if(abs(rts[j])>mx, mx=abs(rts[j])));
  print(CMN[i], " | ", PY, " | N=", nrm, " | ", vector(#rts,j,abs(rts[j])), " | max=", mx);
);
}
print("");
print("=> for every CM orbit the product of the |Y_i| is |N(Y)| >= 1, so max_i |Y_i| >= 1:");
print("   the companion's error term contains a summand of size max|Y_i|^n * n^(-2) >= n^(-2).");
print("   Irrationality needs d_n^2 * |A_n xi - B_n| -> 0, i.e. max|Y_i| < e^(-2) = ", exp(-2), ".");
print("   The ONLY point of X_0(6) with Y = 0 (x = infinity) is the cusp pair 1/2, 1/3 (W = -17).");
print("");
print("(3) magnetic forms whose poles lie only at the cusps 1/2, 1/3 (W = -17):");
{
for(j=1, 8,
  my(Qp, qv, Xs, T, r, Th, Bv, xi, Xi, Phi, rho);
  Qp = (72*x^2+17*x+1)^j;
  qv = vector(poldegree(Qp)+1, t, polcoeff(Qp, t-1));
  Xs = xibasisW(us, F2, qv, 250);
  T = xitomat(Xs, 250);
  r = scanT(T, 60, 150, 250);
  if(r[1]!=1, print("   Q = cusp^", j, "  dim=", #Xs, "  no magnetic element"); next);
  rho = sum(k=0, #Xs-1, r[2][k+1]*x^(k+1))/Qp;
  Phi = F2*subst(numerator(rho),x,us)/subst(denominator(rho),x,us);
  Xi = Dinv(Phi, NQ-2);
  Th = sum(m=1, NQ-2, polcoeff(Xi,m)/m^2*q^m) + O(q^(NQ-1));
  Bv = peel2(Fs*Th, xs, NA, NQ);
  xi = Bv[NA+1]/Av[NA+1];
  print("   Q = cusp^", j, "  dim=", #Xs, "  MAGNETIC, alpha=", r[2]);
  print("        rho = ", rho);
  print("        Apery limit B_n/A_n at n=250 : ", xi, "   rational? ", denominator(xi)<10^20);
  print("        Theta = D^{-2}Xi ; is F*Theta a rational function of x?  F*Theta coefficients in x:");
  print("        ", vector(8,t,Bv[t]));
);
}
quit;
