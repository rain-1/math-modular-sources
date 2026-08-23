default(parisizemax, 8000000000);
default(realprecision, 210);
X  = 0.528749455595766833775958063859865871341153929733052536724383294721129483492053585773843139079492879372388278098384625303363378516146993688539584517920242029624708336;
Y  = 0.373882325600462025027564267690738097743975502372869341714963411208431218497275792601470901538967885214177439871344967969638652927790836356416010161435347320397748845;
X6 = 1.05749891119153366755191612771973174268230785946610507344876658944225896698410717154768627815898575874477655619676925060672675703229398737707916903584048405924941667;
mf32 = mfinit([32,3,-8],0);
eb = mfeigenbasis(mf32);
print("#newforms in S_3^new(32, chi_-8) = ", #eb);
gg = eb[1];
LL = lfunmf(mf32, gg);
L2 = lfun(LL,2);
print("L(g,2) = ", L2);
print();
print("Y - L2/3        = ", Y - L2/3);
print("X - sqrt(2)*L2/3= ", X - sqrt(2)*L2/3);
print("X6 - 2*sqrt(2)*L2/3 = ", X6 - 2*sqrt(2)*L2/3);
print();
/* how many correct digits */
print("log10|Y - L2/3| = ", log(abs(Y-L2/3))/log(10));
print("log10|X - sqrt2 L2/3| = ", log(abs(X-sqrt(2)*L2/3))/log(10));
print("log10|X6 - 2sqrt2 L2/3| = ", log(abs(X6-2*sqrt(2)*L2/3))/log(10));
print();
/* Conrey labels of the quadratic characters mod 8 and mod 32 */
{
foreach([8,32], nn,
  my(G = znstar(nn,1));
  print("modulus ", nn, ":");
  for(n = 1, nn-1, if(gcd(n,nn)==1,
    my(ch = znconreychar(G,n), o = zncharorder(G,ch));
    if(o<=2,
      my(vals = [], ok4=1, ok8=1, okm8=1);
      for(m=1,4*nn, if(gcd(m,nn)==1,
        my(z = chareval(G, ch, m));
        my(v = if(z==0,1,-1));
        if(v != kronecker(-4,m), ok4=0);
        if(v != kronecker(8,m), ok8=0);
        if(v != kronecker(-8,m), okm8=0)));
      print("   Conrey ", nn, ".", n, "  order ", o, "  cond ", zncharconductor(G,ch),
            "  = ", if(o==1,"trivial", if(ok4,"chi_-4", if(ok8,"chi_8", if(okm8,"chi_-8","?")))))))));
}
quit;
