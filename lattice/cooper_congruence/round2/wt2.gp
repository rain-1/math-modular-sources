\\ wt2.gp -- E_4 by reduction, and the weight -2 weakly holomorphic form f = 1/(x F)
\\ on Gamma_0(7), together with its Maass raising  R_{-2} f = -4 pi (Df + f/(2 pi y)).
read("e2.gp");
{ E4hol(z,PR) = my(q=exp(2*Pi*I*z), s=1.); for(n=1,PR, s += 240*sigma(n,3)*q^n); s; }
{ E4(t) = my(R=sl2red(t), z, c, d, PR);
  z=R[1]; c=R[4]; d=R[5];
  PR = ceil(80/(-log(abs(exp(2*Pi*I*z)))/log(10))) + 40;
  E4hol(z,PR)/(c*t+d)^4;
}
{ E2(t) = E2star(t) + 3/(Pi*imag(t)); }
{ uval7(t) = (eta(7*t,1)/eta(t,1))^4; }
{ fhat7(t) = my(u,g,x,F,Dx,DF,f,Df);
  u = uval7(t);
  g = 1 + 13*u + 49*u^2;
  x = u/g;
  F = (7*E2star(7*t) - E2star(t))/6;
  Dx = F*u*(1-49*u^2)/g^2;
  DF = (49*(E2(7*t)^2 - E4(7*t)) - (E2(t)^2 - E4(t)))/72;
  f = 1/(x*F);
  Df = -f^2*(F*Dx + x*DF);
  [f, Df + f/(2*Pi*imag(t))];
}
