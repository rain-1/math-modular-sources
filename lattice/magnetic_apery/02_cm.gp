default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 120);
\\ Apery host: N=6, u = eta1^-5 eta2^1 eta3^-1 eta6^5, C=72
Cc = 72;
uval(t) = eta(t,1)^(-5)*eta(2*t,1)^1*eta(3*t,1)^(-1)*eta(6*t,1)^5;
KMAX = 60;
{
DL = [-3,-4,-7,-8,-11,-12,-15,-16,-19,-20,-23,-24,-27,-28,-32,-35,-36,-39,-40,-43,-48,-51,-52,-55,-56,-64,-67,-72,-75,-84,-88,-96,-100,-120,-147,-163];
for(t=1,#DL,
  my(D); my(vals); my(ims); my(nf);
  D = DL[t];
  vals = List(); ims = List(); nf = 0;
  for(k=1, KMAX,
    my(a); a = 6*k;
    for(b = -a+1, a,
      if((b^2-D) % (4*a) != 0, next);
      my(c); c = (b^2-D)/(4*a);
      if(gcd(gcd(a,b),c) != 1, next);
      nf = nf + 1;
      my(tau); my(uu); my(found);
      tau = (-b + sqrt(D)*I)/(2*a);
      if(D<0, tau = (-b + sqrt(-D)*I)/(2*a));
      uu = uval(tau);
      found = 0;
      for(j=1, #vals,
        if(abs(uu - vals[j]) < 1e-60, found = j; break);
      );
      if(found == 0, listput(vals, uu); listput(ims, imag(tau)),
         if(imag(tau) > ims[found], ims[found] = imag(tau));
      );
    );
  );
  if(#vals == 0, print("D=",D,"  no forms"); next);
  my(pw); my(ws); my(pu);
  ws = vector(#vals, j, vals[j] + 1/(Cc*vals[j]));
  pw = prod(j=1, #ws, (X - ws[j]));
  pu = prod(j=1, #vals, (X - vals[j]));
  print("D=", D, "  #orbits=", #vals, "  maxIm=", vector(#ims,j,ims[j]));
  print("   minpoly_w = ", bestappr(real(pw), 10^40));
  print("   minpoly_u = ", bestappr(real(pu), 10^40));
  print("   w-values  = ", vector(#ws,j,ws[j]));
);
}
quit;
