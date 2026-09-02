\\ 58_s18beta.gp -- s18: try EVERY admissible beta class for each d, and report which give a
\\ real trace (the smallest class is not conjugation-stable, which is why the first run gave
\\ non-integers).
read("56_cdlib.gp");
default(realprecision, 50);
initfser(3,2200);
{
for(d=1,40,
  my(D=-36*d, bs=betas(D,18), row=List());
  if(#bs==0, next);
  for(j=1,#bs,
    my(R=Traw(3,d,bs[j],36), v);
    if(type(R)=="t_STR", listput(row,[bs[j],"NOREP"]); next);
    v = R[1]/sqrt(3);
    listput(row, [bs[j],
      if(abs(imag(v))>1e-20*(1+abs(v)), "CPLX",
         if(abs(real(v)-round(real(v)))<1e-20, round(real(v)), "NONINT"))]));
  print("d=",d,"  betas ",bs,"   -> ",Vec(row));
);
}
quit;
