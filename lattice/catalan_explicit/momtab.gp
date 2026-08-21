\p 400
G=Catalan;
{ for(n=1,5, my(m=3*n, D=lcm(vector(6*n,i,i)), ee2=ee(m));
  print("\n=== n=",n,"  m=3n=",m,"   log D_{6n}^2 /n = ",log(1.0*D^2)/n,"  e_{3n}=",ee2," ===");
  print(" j | v2(denA) v2(denB) | log10 denA log10 denB | log|M|/n | log(den*|M|)/n | log2(den) ");
  for(j=0,m, my(r=nestgen(m,j), A=r[1], B=r[2], dA=denominator(A), dB=denominator(B), d=lcm(dA,dB));
    my(val=A*G-B);
    printf("%3d | %7d %7d | %10.3f %10.3f | %9.4f | %9.4f | %9.3f\n",
      j, valuation(dA,2), valuation(dB,2), log(1.0*dA)/log(10), log(1.0*dB)/log(10),
      log(abs(val))/n, log(1.0*d*abs(val))/n, log(1.0*d)/log(2));
  );
); }
\q
