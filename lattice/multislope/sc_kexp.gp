default(parisizemax, 6000000000);
read("sc_rows.gp");
{
rows = [["R1",R1cf,2,1500],["R2",R2cf,2,1500],["R3",R3cf,2,1500],["R4",R4cf,6,1000],["R5",R5cf,4,1000]];
for(i=1,#rows,
  my(nm=rows[i][1], cf=rows[i][2], r=rows[i][3], N=rows[i][4]);
  my(ks = vector(r-1, j, denexp(compan(cf,r,j,N), N)));
  print(nm, ": N=", N, "  denominator exponent per companion = ", ks);
);
}
