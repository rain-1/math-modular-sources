default(parisizemax, 6000000000);
read("sc_rows.gp");
{
A4 = r4A(200);
for(ord=2,5,
  for(d=2,14,
    my(f = scfit(A4, ord, d));
    if(type(f)!="t_INT" || f!=0, print("ord ",ord," d ",d," -> ", if(type(f)=="t_VEC","VEC",f)));
  );
);
print("scan done");
}
