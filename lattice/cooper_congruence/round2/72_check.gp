\\ 72_check.gp -- independent checks:
\\ (1) S_4(Gamma_0(7)) and its Fricke eigenvalue (the uniqueness input at level 28);
\\ (2) the NEW divisibility  p^{2a} | c(p^{2a} d)  on the extended c(d) table, INCLUDING
\\     the cells with p | d, which the PZ agent could not test.
print("=== (1) S_4(Gamma_0(7)) ===");
S4 = mfinit([7,4],1);
print("dim S_4(Gamma_0(7)) = ", mfdim(S4));
print("dim M_4(Gamma_0(7)) = ", mfdim(mfinit([7,4],4)));
{ my(F=mfeigenbasis(S4), v);
  v = mfcoefs(F[1], 15);
  print("  eigenform coefficients q^0..q^15: ", v);
  print("  a_7 = ", v[8]);
  print("  Atkin-Lehner W_7 eigenvalue on S_4: ", mfatkineigenvalues(S4,7)); }
print("  => dim S_4(Gamma_0(7))^{W_7=-1} = 0   (Phi|W_7 = -Phi lives there)");
print("");
print("=== (2) p^{2a} | c(p^{2a} d) on the d <= 1600 table ===");
L = read("/home/ubuntu/code/math-modular-sources/lattice/cooper_congruence/round2/71_cd_s7_big.txt");
print("table entries: ", #L, "   max d: ", L[#L][1]);
CD = Map();
{ for(i=1,#L, if(type(L[i][2])=="t_INT", mapput(CD, L[i][1], L[i][2]))); }
{ my(tot=0, bad=0);
  forprime(p=2,23,
    for(a=1,4,
      my(P=p^(2*a), cn=0, cnp=0, mn=99, mnp=99);
      for(d=1,1600,
        if(!mapisdefined(CD,d) || !mapisdefined(CD,P*d), next);
        my(v=mapget(CD,P*d));
        if(v==0, next);
        tot++;
        if(v % P != 0, bad++; print("   *** FAIL p=",p," a=",a," d=",d,"  c=",v));
        my(e=valuation(v,p)-2*a);
        if(d%p==0, cnp++; if(e<mnp,mnp=e), cn++; if(e<mn,mn=e)));
      if(cn+cnp>0,
        print("  p=",p," a=",a,":  p nmid d: ",cn," tests, min excess ",if(cn>0,mn,"-"),
              "   |   p | d: ",cnp," tests, min excess ",if(cnp>0,mnp,"-"))));
  );
  print("  TOTAL ",tot," tests, ",bad," failures");
}
quit;
