default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 150);
Cc = 72;
uval(t) = eta(t,1)^(-5)*eta(2*t,1)^1*eta(3*t,1)^(-1)*eta(6*t,1)^5;
KMAX = 80;
prim(P) = my(d,g); d=1; for(j=0,poldegree(P), d=lcm(d,denominator(polcoeff(P,j)))); P=P*d; g=content(P); P/g;
{
DL = [-8,-12,-15,-20,-23,-24,-32,-35,-36,-39,-44,-47,-48,-56,-59,-71,-72,-80,-83,-95,-96,-104,-107,-119,-120];
print("CM points on X_0(6):  D | #orbits | maxIm | primitive integer minpoly of w=u+1/(72u) | is it a square?");
for(t=1,#DL,
  my(D, vals, ims, ws, pw, PW, mi, sq);
  D = DL[t];
  vals = List(); ims = List();
  for(k=1, KMAX,
    my(a); a = 6*k;
    for(b = -a+1, a,
      if((b^2-D) % (4*a) != 0, next);
      my(c, tau, uu, found);
      c = (b^2-D)/(4*a);
      if(gcd(gcd(a,b),c) != 1, next);
      tau = (-b + sqrt(-D)*I)/(2*a);
      uu = uval(tau);
      found = 0;
      for(j=1, #vals, if(abs(uu - vals[j]) < 1e-80, found = j; break));
      if(found == 0, listput(vals, uu); listput(ims, imag(tau)), if(imag(tau) > ims[found], ims[found] = imag(tau)));
    );
  );
  if(#vals == 0, print("D=",D,"  no forms"); next);
  ws = List();
  for(j=1,#vals, my(wv, fnd); wv = vals[j] + 1/(Cc*vals[j]); fnd=0; for(t2=1,#ws, if(abs(wv-ws[t2])<1e-80, fnd=1; break)); if(fnd==0, listput(ws, wv)));
  pw = prod(j=1, #ws, (X - ws[j]));
  PW = prim(bestappr(real(pw), 10^50));
  mi = 0; for(j=1,#ims, if(ims[j]>mi, mi=ims[j]));
  sq = issquare(PW);
  print("D=", D, " | ", #vals, " orbits, ", #ws, " w-values | maxIm=", mi, " | ", PW, " | sq=", sq);
);
}
quit;
