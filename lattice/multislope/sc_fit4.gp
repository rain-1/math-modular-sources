default(parisizemax, 6000000000);
read("sc_rows.gp");
{
A4 = r4A(150);
print("A_0..A_8 = ", vector(9,i,A4[i]));
found = 0;
for(ord=2,6,
  for(d=2,8,
    if(found, next);
    my(f = scfit(A4, ord, d));
    if(type(f)=="t_VEC",
      found = 1;
      print("order ", ord, " deg ", d, "  FOUND");
      for(i=1,#f, print("  P_",i-1," = ", f[i]));
      my(bad=0);
      for(nn=ord, 148, my(s=0); for(i=0,ord, s += subst(f[i+1],'n,nn)*A4[nn-i+1]); if(s!=0, bad++));
      print("  verification failures on n=",ord,"..148: ", bad);
      my(dg = poldegree(f[1]));
      my(lead = vector(#f, i, polcoeff(f[i], dg, 'n)));
      print("  leading coeffs: ", lead);
      my(cp = sum(i=1,#f, lead[i]*'x^(#f-i)));
      print("  char poly = ", factor(cp));
    );
  );
);
}
