\\ 66_thetalift.gp -- the full theta-lift prediction at level one:
\\   a(d) = sqrt(d) * Tr_{-3d}(fhat) / 192   for EVERY d = 0,1 mod 4,
\\ where a(d) are the coefficients of Pasol-Zudilin's weight-5/2 form
\\   f_4a = (7/8) g0 + (1/768) g1 - (1/768) g2 = q^{-3} + q - 506 q^4 + ...
read("heeg.gp"); read("e2.gp");
default(realprecision, 50);
M = 200;
q4 = 'q^4;
{ Es(k,d,MM) = my(s); s = 1 + O('q^(MM+1));
  if(k==2, s = 1 - 24*sum(n=1,MM\d, sigma(n)*'q^(d*n)));
  if(k==4, s = 1 + 240*sum(n=1,MM\d, sigma(n,3)*'q^(d*n)));
  if(k==6, s = 1 - 504*sum(n=1,MM\d, sigma(n,5)*'q^(d*n)));
  s; }
th  = 1 + 2*sum(n=1,sqrtint(M),'q^(n^2)) + O('q^(M+1));
E24 = (-Es(2,1,M) + 3*Es(2,2,M) - 2*Es(2,4,M))/24;
D4  = 'q^4*prod(n=1,M\4,(1-'q^(4*n)+O('q^(M+1)))^24);
g0 = th*(th^4 - 20*E24);
g1 = th*Es(4,4,M)^2*Es(6,4,M)/D4;
g2 = g0*Es(4,4,M)^3/D4;
f4a = (7/8)*g0 + (1/768)*g1 - (1/768)*g2;
print("f_4a coefficients q^-3 .. q^13: ", vector(17,i,polcoeff(f4a,i-4)));
\\ fhat on SL_2(Z)
E4s = Es(4,1,M); E6s = Es(6,1,M); Dl = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
fw = E4s*E6s/Dl; Dfw = 'q*deriv(fw,'q);
{ evs(f,q,TT) = my(v=valuation(f,'q), s=0.); for(n=v,TT, s += polcoeff(f,n)*q^n); s; }
{ fhat1(t) = my(R=sl2red(t), z, qq, y); z=R[1]; y=imag(z); qq=exp(2*Pi*I*z);
  evs(Dfw,qq,M-3) + evs(fw,qq,M-3)/(2*Pi*y); }
{ trD(D) = my(RF=redforms(D), t=0.);
  for(i=1,#RF, t += genchar(RF[i],-3)*fhat1((-RF[i][2]+I*sqrt(-D))/(2*RF[i][1]))/omeg(RF[i]));
  t; }
print("");
print(" d   a(d)          sqrt(d)*Tr_{-3d}(fhat)/192      difference");
{
for(d=1,40,
  if(d%4!=0 && d%4!=1, next);
  my(a=polcoeff(f4a,d), T=trD(-3*d), P=sqrt(d)*T/192);
  print("  ",d,"   ",a,"    ",P,"    ",P-a);
);
}
quit;
