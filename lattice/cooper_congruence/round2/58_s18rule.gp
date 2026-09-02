\\ 58_s18rule.gp -- s18 with the SQUARE-RULE beta (18 for d odd, 0 for d even, which is
\\ 18m mod 36 at d = m^2): classify T(d)/sqrt3 as ZERO / INTEGER / NOT-RATIONAL, and
\\ correlate with d mod 3 and with d being a square.
read("56_cdlib.gp");
default(realprecision, 50);
initfser(3,2200);
Z=List(); INT=List(); IRR=List();
{
for(d=1,120,
  my(bt=if(d%2, 18, 0), D=-36*d, R, v);
  if(!setsearch(Set(betas(D,18)), bt), print("  d=",d," square-rule beta ",bt," NOT ADMISSIBLE"); next);
  R = Traw(3,d,bt,36);
  if(type(R)=="t_STR", print("  d=",d," NOREP"); next);
  v = R[1]/sqrt(3);
  if(abs(v) < 1e-20, listput(Z,d),
     if(abs(imag(v))<1e-18*(1+abs(v)) && abs(real(v)-round(real(v)))<1e-18, listput(INT,[d,round(real(v))]),
        listput(IRR,d)));
);
}
print("ZERO   (",#Vec(Z)," values): ",Vec(Z));
print("   all = 2 mod 3 ? ", vecmin(vector(#Vec(Z),i,Vec(Z)[i]%3==2))==1);
print();
print("INTEGER nonzero (",#Vec(INT)," values): ",Vec(INT));
print("   all squares ? ", vecmin(vector(#Vec(INT),i,issquare(Vec(INT)[i][1])))==1);
print();
print("NOT RATIONAL (",#Vec(IRR)," values): ",Vec(IRR));
print("   any square among them ? ", if(#Vec(IRR)>0, vecmax(vector(#Vec(IRR),i,issquare(Vec(IRR)[i]))), 0));
print("   d mod 3 values present: ", Set(vector(#Vec(IRR),i,Vec(IRR)[i]%3)));
quit;
