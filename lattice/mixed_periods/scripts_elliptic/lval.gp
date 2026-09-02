default(parisize,"4G");
default(realprecision,100);
doit(name,av)=
{
  my(p,a,b,c,d,E,Em,N);
  p = prod(i=1,#av,(1-av[i]*x));
  a=polcoeff(p,3); b=polcoeff(p,2); c=polcoeff(p,1);
  E=ellinit([0,b,0,a*c,a^2*1]); Em=ellminimalmodel(E); N=ellglobalred(E)[1];
  print(name,"|",N,"|",lfun(Em,0,1),"|",lfun(Em,2),"|",Em.omega[1],"|",imag(Em.omega[2]),"|",E.omega[1],"|",imag(E.omega[2]),"|",ellanalyticrank(Em)[1],"|",E.j);
}
for(k=1,7, my(mv=[1,2,3,4,5,6,12], m=mv[k]); doit(Str("E_",m), [1,9,4*m]));
for(m=1,3, doit(Str("F_",m), [1,25,4*m]));
doit("G",[4,8,12]);
quit;
