\\ 57_cd_s10.gp -- c(d) for row s10 over ALL admissible d.  D = -4d, chi_{-4},
\\ beta = 6*sqrt(d) mod 20 on squares, else the smallest admissible class.
\\ Raw output: A(d) := T(d)/i  (T = sum chi_{-4}(Q) fhat/omega).  Normalisations are
\\ applied afterwards:  c(d) = i*A(d)*(-i) = A(d)  with lambda = i, or 2*A(d) with lambda = 2i.
read("56_cdlib.gp");
default(realprecision, 60);
initfser(2,2200);
DMAX = 400;
L = List();
gettime();
{
for(d=1,DMAX,
  if(d%4!=0 && d%4!=1, next);            \\ complementary discriminant D1 = d
  my(bt=chosenbeta(2,d), R, v);
  if(bt<0, next);
  R = Traw(2,d,bt,36);
  if(type(R)=="t_STR", listput(L,[d,bt,0,"NOREP"]); next);
  v = R[1]/I;
  if(abs(imag(v)) > 1e-25*(1+abs(v)),
     listput(L,[d,bt,R[2],"DEGEN"]),
     listput(L,[d,bt,R[2],round(real(v))*if(abs(real(v)-round(real(v)))<1e-25,1,"NONINT")]));
);
}
print("s10: ", #Vec(L), " admissible d <= ", DMAX, "   [", gettime(), " ms]");
write("57_cd_s10.txt", Vec(L));
print(vector(min(30,#Vec(L)), i, Vec(L)[i]));
quit;
