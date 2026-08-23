default(parisizemax, 2000000000);
{
foreach([8,32], nn,
  my(G = znstar(nn,1));
  print("modulus ", nn, ":");
  for(n = 1, nn-1, if(gcd(n,nn)==1,
    my(ch = znconreylog(G,n), o = charorder(G,ch));
    if(o<=2,
      my(ok4=1, ok8=1, okm8=1);
      for(m=1,8*nn, if(gcd(m,nn)==1,
        my(z = chareval(G, ch, m), v);
        v = if(z==0,1,-1);
        if(v != kronecker(-4,m), ok4=0);
        if(v != kronecker(8,m), ok8=0);
        if(v != kronecker(-8,m), okm8=0)));
      print("   Conrey ", nn, ".", n, "  order ", o, "  cond ", znconreyconductor(G,ch),
            "  = ", if(o==1,"trivial", if(ok4,"chi_-4", if(ok8,"chi_8", if(okm8,"chi_-8","?")))))))));
}
quit;
