\\ 13_ident7.gp -- THE IDENTIFICATION (row s7):
\\   beta(m) = i sqrt(3) * sum_Q chi_{-3}(Q) omega_Q^{-1} (R_{-2} f / (-4 pi))(alpha_Q)
\\ over the Heegner points of discriminant -3 m^2 on Gamma_0(7) with B = 5m mod 14,
\\ where f = 1/(x F) is THE weakly holomorphic weight -2 form on Gamma_0(7) with
\\ q-expansion q^{-1} + O(1) and f|W_7 = -f, and R_{-2} f = -4 pi (Df + f/(2 pi y)).
read("lib.gp"); read("heeg.gp"); read("wt2b.gp");
default(realprecision, 60);
N=7; D0=-3;
{ tr(m) = my(d=-3*m^2, bt, RF, rep, al, ch, om, t=0., nf=0);
  bt = (5*m)%14;
  RF = redforms(d);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, print("NOREP ",m," ",RF[i]); next);
    ch = genchar(rep,D0); om = omeg(rep); nf++;
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    t += ch*fhatq(al,7,FS7,400)/om);
  [t,nf,#RF];
}
bet = read("beta_s7.txt");
cp  = read("cp_s7.txt");
{
for(m=1,26,
  my(T=tr(m), P);
  P = I*sqrt(3)*T[1];
  print("m=",m,"  #cls=",T[3],"/",T[2],"   i sqrt3 * Tr = ", P, "     beta = ", bet[m], "     diff = ", P-bet[m]);
);
}
quit;
