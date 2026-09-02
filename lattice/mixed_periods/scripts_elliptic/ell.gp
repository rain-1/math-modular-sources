default(parisize,"2G");
default(realprecision,60);
doit(name,av)=
{
  my(p,a,b,c,d,E,Em,N,j,per,L0,L0p,L2,Lstar);
  p = prod(i=1,#av,(1-av[i]*x));
  a = polcoeff(p,3); b = polcoeff(p,2); c = polcoeff(p,1); d = polcoeff(p,0);
  \\ y^2 = a x^3 + b x^2 + c x + d ;  x=X/a, y=Y/a  => Y^2 = X^3 + b X^2 + a c X + a^2 d
  E = ellinit([0,b,0,a*c,a^2*d]);
  Em = ellminimalmodel(E);
  N = ellglobalred(E)[1];
  j = E.j;
  per = Em.omega;
  print("### ", name, "   P = ", p);
  print("  Weierstrass (from x=X/a,y=Y/a): [0,",b,",0,",a*c,",",a^2*d,"]");
  print("  minimal model: ", Em.a1,",",Em.a2,",",Em.a3,",",Em.a4,",",Em.a6, "   disc=",Em.disc);
  print("  conductor N = ", N, "   j = ", j, "   (j factored: ", if(type(j)=="t_FRAC"||type(j)=="t_INT", factor(numerator(j))," "), ")");
  print("  ellrank/analytic rank: ", ellanalyticrank(Em)[1]);
  print("  omega (minimal model) = ", Em.omega);
  print("  omega (our model)     = ", E.omega);
  L2 = lfun(Em,2);
  L0p = lfun(Em,0,1);
  L0 = lfun(Em,0);
  print("  L(E,2)   = ", L2);
  print("  L(E,0)   = ", L0);
  print("  L'(E,0)  = ", L0p);
  print("  L(E,1)   = ", lfun(Em,1));
  print("  L'(E,1)  = ", lfun(Em,1,1));
  print("");
}
for(k=1,7, my(mv=[1,2,3,4,5,6,12], m=mv[k]); doit(Str("E_",m), [1,9,4*m]));
for(m=1,3, doit(Str("F_",m), [1,25,4*m]));
doit("G", [4,8,12]);
quit;
