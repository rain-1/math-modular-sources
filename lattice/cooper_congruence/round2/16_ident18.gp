\\ 16_ident18.gp -- identification attempt for row s18 (N=18, pole disc -36).
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 50);
N=18; initfser(3,700);
{ tr(m, b0, D0) = my(d=-36*m^2, bt, RF, rep, al, ch, om, t=0., nf=0);
  bt = (b0*m)%(2*N);
  RF = redforms(d);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,80);
    if(rep==0, next);
    ch = if(D0==0, 1, genchar(rep,D0)); om = omeg(rep); nf++;
    al = (-rep[2] + I*m*6)/(2*rep[1]);
    t += ch*fhatR(3,al)/om);
  [t,nf,#RF];
}
bet = read("beta_s18.txt");
{
for(m=1,10,
  my(A=tr(m,18,-3), B=tr(m,18,-4), C=tr(m,18,0));
  print("m=",m," cls ",A[3],"/",A[2],"   beta=",bet[m]);
  print("    chi_-3: ", A[1], "   chi_-4: ", B[1], "   nochi: ", C[1]);
);
}
quit;
