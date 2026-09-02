\\ e2.gp -- E_2^*(tau) = E_2(tau) - 3/(pi*Im tau), weight 2 for SL_2(Z), by reduction.
{ sl2red(t) = my(a=1,b=0,c=0,d=1, n, z=t, ab);
  for(it=1,200,
    n = round(real(z));
    z = z - n; a -= n*c; b -= n*d;   \\ z <- z-n ;  matrix update for gamma with gamma*t = z
    ab = abs(z);
    if(ab >= 1-1e-30, break);
    z = -1/z;
    [a,b,c,d] = [c,d,-a,-b];
  );
  [z,a,b,c,d];
}
{ E2hol(z, PR) = my(q=exp(2*Pi*I*z), s=1., t);
  for(n=1, PR, s -= 24*sigma(n)*q^n);
  s;
}
{ E2star(t) = my(R=sl2red(t), z, c, d, v);
  z = R[1]; c = R[4]; d = R[5];
  v = E2hol(z, ceil(60/(-log(abs(exp(2*Pi*I*z))) /log(10))) + 40 );
  v = v - 3/(Pi*imag(z));
  \\ E2*(t) = (c t + d)^{-2} E2*(gamma t) with gamma t = z
  v/(c*t+d)^2;
}
