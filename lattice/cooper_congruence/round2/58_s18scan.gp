\\ 58_s18scan.gp -- s18: all beta classes, d = 41..150, showing the RAW value B = T/sqrt(3)
\\ so that the "NONINT" entries can be identified.
read("56_cdlib.gp");
default(realprecision, 40);
initfser(3,2200);
{
for(d=41,150,
  my(D=-36*d, bs=betas(D,18), row=List(), any=0);
  if(#bs==0, next);
  for(j=1,#bs,
    my(R=Traw(3,d,bs[j],36), v);
    if(type(R)=="t_STR", listput(row,[bs[j],"NOREP"]); next);
    v = R[1]/sqrt(3);
    if(abs(v)>1e-15, any=1);
    listput(row, [bs[j], v]));
  if(!any, print("d=",d,"  ALL ZERO   ",if(issquare(d),"(SQUARE)","")); next);
  print("d=",d,if(issquare(d),"  (SQUARE)",""),"  betas ",bs);
  for(j=1,#row, print("     b=",Vec(row)[j][1],"  B=",Vec(row)[j][2]));
);
}
quit;
