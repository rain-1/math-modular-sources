\\ 03_fit7.gp -- trace sequences on X_0(7) at disc -3m^2, one beta-class, all forms,
\\ for the basis {1, u, 1/u} of functions with poles of order <=1 at the two cusps.
read("lib.gp"); read("heeg.gp");
default(realprecision, 120);
N = 7;
{ tr7(m, PRIM) = my(d=-3*m^2, bt, RF, rep, al, om, t1, tu, ti, g, w, nf);
  bt = (5*m)%14;
  RF = redforms(d);
  t1=0.; tu=0.; ti=0.; nf=0;
  for(i=1,#RF,
    g = gcd(gcd(RF[i][1],RF[i][2]),RF[i][3]);
    if(PRIM && g>1, next);
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, print("  NOREP m=",m," ",RF[i]); next);
    nf++;
    om = omeg(rep);
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    w = (eta(7*al,1)/eta(al,1))^4;
    t1 += 1./om; tu += w/om; ti += (1/w)/om);
  [t1,tu,ti,nf];
}
bet = read("beta_s7.txt");
{
print("m  [T1, Tu, T1/u]  (all forms), then (primitive only)");
for(m=1,20,
  my(A=tr7(m,0), Bp=tr7(m,1));
  print(m, "  all: ", A[4], " ", A[1], " | ", A[2], " | ", A[3]);
  print("   prim: ", Bp[4], " ", Bp[1], " | ", Bp[2], " | ", Bp[3]);
  print("   beta=",bet[m]);
);
}
quit;
