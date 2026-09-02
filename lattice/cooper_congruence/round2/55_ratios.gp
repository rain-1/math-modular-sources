\\ 55_ratios.gp -- consolidated ratio test for the three rows.
\\ Claim under test:  beta(m)/tau(m) = L*(1 - kappa/m) + O(R^{-m/2}),  kappa = 1/(2 pi Im tau_0).
\\ If beta were a linear combination of twisted CM traces of modular FUNCTIONS,
\\ the ratio would be L + O(R^{-m/2}) with no 1/m term.
read("50_lib.gp");
default(realprecision, 60);
BET = [read("20_beta_s7.txt"), read("20_beta_s10.txt"), read("20_beta_s18.txt")];
TAU = [read("52_tau_s7.txt"), read("53_tau_s10.txt"), read("54_tau_s18.txt")];
KAP = [7/(Pi*sqrt(3)), 5/Pi, 3/Pi];
RRV = [exp(Pi*sqrt(3)/7), exp(Pi/5), exp(Pi/3)];
\\ empirical limit L, by class of m
{ Lof(k,m) =
  if(k==1, return(3));
  if(k==2, return(if(m%2==1, 2, 4)));
  if(m%2==0, return(0));
  if(m%3==0, 1/2, 1/4);
}
{
for(k=1,3,
  my(bet=BET[k], tt=TAU[k], kap=KAP[k], R=RRV[k], MM=#tt);
  print("======== row ", NAM[k], "   kappa=",kap,"   R=",R);
  print("   m     L    beta/tau                      r1-L (should ~ -L kappa/m)      r2=m beta/((m-kappa)tau)         (r2-L)*R^(m/2)");
  for(m=1,MM,
    if(tt[m]==0, print("   ",m,"    anomalous / tau=0"); next);
    my(L=Lof(k,m), r1, r2);
    if(L==0, print("   ",m,"    anomalous"); next);
    r1 = bet[m]/tt[m]*1.;
    r2 = m*bet[m]/((m-kap)*tt[m]);
    print("   ",m,"   ",L,"   ",r1,"   ",(r1-L)*m," vs ",-L*kap,"   ",r2,"   ",(r2-L)*R^(m/2)));
  print();
);
}
print();
print("== SUMMARY: sup over the last 10 clean m of |(r2-L)*R^(m/2)|  and  of |(r1-L)*R^(m/2)| ==");
{
for(k=1,3,
  my(bet=BET[k], tt=TAU[k], kap=KAP[k], R=RRV[k], MM=#tt, s2=0., s1=0., c=0);
  forstep(m=MM,1,-1,
    if(c>=10, break);
    if(tt[m]==0, next);
    my(L=Lof(k,m)); if(L==0, next);
    c++;
    s2 = max(s2, abs(m*bet[m]/((m-kap)*tt[m]) - L)*R^(m/2));
    s1 = max(s1, abs(bet[m]/tt[m]*1. - L)*R^(m/2)));
  print("row ",NAM[k],":  with the (m-kappa) correction: ",s2,"     WITHOUT it: ",s1);
);
}
quit;
