\\ 74_hecke.gp -- on the extended table (d <= 2500):
\\ (A) the divisibility p^{2a} | c(p^{2a} d), including the p | d cells;
\\ (B) the Hecke fit  c(p^2 d) + tau*p*c(d) + p^3 c(d/p^2) = lambda c(d) mod p^2,
\\     separately for p nmid d and p | d, and for tau = 1 versus tau = 0 (the value a
\\     genuine T_{p^2} would take when p | d).
L = read("/home/ubuntu/code/math-modular-sources/lattice/cooper_congruence/round2/73_cd_s7_2500.txt");
CD = Map();
{ for(i=1,#L, if(type(L[i][2])=="t_INT", mapput(CD, L[i][1], L[i][2]))); }
get(d) = if(mapisdefined(CD,d), mapget(CD,d), "NA");
print("table entries ", #L, "  max d ", L[#L][1]);
print("");
print("=== (A)  p^{2a} | c(p^{2a} d) ===");
{ my(tot=0, bad=0);
  forprime(p=2,29,
    for(a=1,4,
      my(P=p^(2*a), cn=0, cnp=0, mn=99, mnp=99);
      for(d=1,2500,
        if(get(d)=="NA" || get(P*d)=="NA", next);
        my(v=get(P*d)); if(v==0, next);
        tot++;
        if(v%P!=0, bad++; print("  *** FAIL p=",p," a=",a," d=",d));
        my(e=valuation(v,p)-2*a);
        if(d%p==0, cnp++; if(e<mnp,mnp=e), cn++; if(e<mn,mn=e)));
      if(cn+cnp>0, print("  p=",p," a=",a,":  p nmid d ",cn," tests, min excess ",if(cn>0,mn,"-"),
                          "  |  p | d ",cnp," tests, min excess ",if(cnp>0,mnp,"-")))));
  print("  TOTAL ",tot," tests, ",bad," failures");
}
print("");
print("=== (B)  lambda mod p^2 from  c(p^2 d) + tau p c(d) + p^3 c(d/p^2) = lambda c(d) ===");
{ forprime(p=2,23,
    my(S1=Set(), S0=Set(), T1=Set(), T0=Set(), n1=0, n0=0);
    for(d=1,2500,
      if(get(d)=="NA" || get(p^2*d)=="NA", next);
      my(cd=get(d), c2=get(p^2*d), cdp = if(d%(p^2)==0 && get(d/p^2)!="NA", get(d/p^2), 0), l1, l0);
      if(cd==0 || cd%p==0, next);
      l1 = lift(Mod(c2 + p*cd + p^3*cdp, p^2)/Mod(cd,p^2));
      l0 = lift(Mod(c2 + p^3*cdp, p^2)/Mod(cd,p^2));
      if(d%p==0, n0++; T1=setunion(T1,Set([l1])); T0=setunion(T0,Set([l0])),
                 n1++; S1=setunion(S1,Set([l1])); S0=setunion(S0,Set([l0]))));
    print("  p=",p,"  p should be ",p);
    print("     p nmid d (",n1," tests):  tau=1 -> ",S1,"   tau=0 -> ",S0);
    print("     p |  d   (",n0," tests):  tau=1 -> ",T1,"   tau=0 -> ",T0));
}
quit;
