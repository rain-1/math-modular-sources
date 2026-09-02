\\ 21_odd.gp -- periodicity of gamma on the odd-index subsequence (mod 4, 8, 16, 32).
read("lib.gp");
G  = [read("20b_gamma_s7.txt"), read("20b_gamma_s10.txt"), read("20b_gamma_s18.txt")];
N  = #G[1];
print("N = ", N);
{
for(k=1,3,
  my(W = vector(N\2, i, G[k][2*i-1]));
  for(e=2,5,
    my(m=2^e, T=0, s);
    for(t=1,600,
      my(ok=1);
      for(i=1,#W-t, if((W[i]-W[i+t])%m!=0, ok=0; break));
      if(ok, T=t; break);
    );
    s = if(T>0, concat("period ",T), "no period <= 600");
    print("  ",NAM[k],"  gamma(2i-1) mod ",m,":  ",s);
  );
);
}
quit;
