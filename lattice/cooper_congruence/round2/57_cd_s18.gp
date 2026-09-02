\\ 57_cd_s18.gp -- c(d) for row s18 over ALL admissible d.  D = -36d, chi* = chi_{-3}*(-4|cont),
\\ admissible forms 3 nmid cont, beta = 18*sqrt(d) mod 36 on squares else smallest class.
\\ Raw output: B(d) := T(d)/sqrt(3).
read("56_cdlib.gp");
default(realprecision, 60);
initfser(3,2200);
DMAX = 300;
L = List();
gettime();
{
for(d=1,DMAX,
  my(bt=chosenbeta(3,d), R, v);
  if(bt<0, next);
  R = Traw(3,d,bt,36);
  if(type(R)=="t_STR", listput(L,[d,bt,0,"NOREP"]); next);
  v = R[1]/sqrt(3);
  if(abs(imag(v)) > 1e-22*(1+abs(v)),
     listput(L,[d,bt,R[2],"DEGEN"]),
     listput(L,[d,bt,R[2],if(abs(real(v)-round(real(v)))<1e-22, round(real(v)), "NONINT")]));
);
}
print("s18: ", #Vec(L), " admissible d <= ", DMAX, "   [", gettime(), " ms]");
write("57_cd_s18.txt", Vec(L));
print(vector(min(30,#Vec(L)), i, Vec(L)[i]));
quit;
