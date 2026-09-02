\\ wt2b.gp -- f = 1/(xF), weight -2 weakly holomorphic on Gamma_0(7), by q-expansion,
\\ with reduction to the fundamental domain of the Fricke group Gamma_0(7)+ (f|W_7 = -f).
read("lib.gp");
NQ = 500;
{ mkf(k) = my(S=Setup(k,NQ), x=S[3], F=S[2], f);
  f = 1/(x*F);
  [f, 'q*deriv(f,'q)];
}
FS7 = mkf(1);
\\ reduce tau by Gamma_0(N)+ ; returns [tau_red, sign]
{ redN(t, N) = my(s=1, n);
  for(it=1,300,
    n = round(real(t)); t = t - n;
    if(N*abs(t)^2 >= 1-1e-40, break);
    t = -1/(N*t); s = -s);
  [t,s];
}
{ evser(f, q, TT) = my(v=valuation(f,'q), s=0.);
  for(n=v, TT, s += polcoeff(f,n)*q^n);
  s;
}
{ fhatq(t, N, FF, TT) = my(R=redN(t,N), z, q, y);
  z = R[1]; y = imag(z); q = exp(2*Pi*I*z);
  R[2]*( evser(FF[2],q,TT) + evser(FF[1],q,TT)/(2*Pi*y) );
}
