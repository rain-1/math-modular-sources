\\ 65_trace7disc.gp -- the trace formula at another discriminant: Psi(f_7), D_0=-7,
\\ level one.  Which character in the trace corresponds to which psi in the Mobius
\\ inversion?  (Expected from the F_4a data: untwisted trace <-> psi = chi_{D_0}.)
read("heeg.gp"); read("e2.gp");
default(realprecision, 50);
M = 200;
E4s = 1 + 240*sum(n=1,M,sigma(n,3)*'q^n) + O('q^(M+1));
E6s = 1 - 504*sum(n=1,M,sigma(n,5)*'q^n) + O('q^(M+1));
Dl  = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
js  = E4s^3/Dl;
fw  = E4s*E6s/Dl;
Dfw = 'q*deriv(fw,'q);
Ph  = 27*E4s*(19*js-8*3375)/(js+3375)^2;
{ betav(D) = my(c=vector(M,m,polcoeff(Ph,m)), cp, b=vector(M));
  cp = vector(M,m,c[m]/m);
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*if(D==0,1,kronecker(D,d))*cp[n/d]); b[n]=s); b; }
B1 = betav(0); Bx = betav(-7);
{ evs(f,q,TT) = my(v=valuation(f,'q), s=0.); for(n=v,TT, s += polcoeff(f,n)*q^n); s; }
{ fhat1(t) = my(R=sl2red(t), z, q, y);
  z=R[1]; y=imag(z); q=exp(2*Pi*I*z);
  evs(Dfw,q,M-3) + evs(fw,q,M-3)/(2*Pi*y); }
{ tr(m,tw) = my(d=-7*m^2, RF=redforms(d), t=0.);
  for(i=1,#RF, t += if(tw,genchar(RF[i],-7),1)*fhat1((-RF[i][2]+I*m*sqrt(7))/(2*RF[i][1]))/omeg(RF[i]));
  t; }
{ my(l1, lx);
  l1 = B1[1]/tr(1,1); lx = Bx[1]/tr(1,0);
  print("lambda (twisted trace  <-> psi=1     ) = ", l1);
  print("lambda (untwisted trace <-> psi=chi_-7) = ", lx);
  for(m=1,12,
    print("m=",m,"   psi=1:   err = ", l1*tr(m,1)-B1[m],
                "     psi=chi_-7:  err = ", lx*tr(m,0)-Bx[m]));
}
quit;
