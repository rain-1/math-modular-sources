default(realprecision, 120);
LK(D, s) = { my(q=abs(D), t=0); for(j=1,q-1, t += kronecker(D,j)*zetahurwitz(s, j/q)); t/q^s; }

xi1 = 0.06786721954287826166646998964563721793157057974785253811627679599873293240387754600472436904131064158841277383301338254691816105592800291978569842592232280720350666634075477524078495373724845256233523035364407420441444991361676976050336957323807362044 - 0.06907819963322225514670799790743491919071785808412459538740013509317488539578332659339145894820192577647267658296914645712671900760850497706419803607235071864675633771938344206632717540610599368109439301877396417892631629513247150535996752102799884096*I;
xi2 = 0.1642087341590574584714225970206947065830122999536287405694316515422586105676807745518762686902555647720236146720735800852953327485050302342736826060045161009106544716484284489377941968395761197077123295366087484842016287741923684399374479825889396610 - 0.2369659751659112010711143299806665676669665318863734518080024379160442279419310559148346593990524891514289658806059918517094079894759328173258155337192594168909099387369602162580142247035701450032542277433112273060653991607390105577213119189860336177*I;

\\ ---- level 7 weight 3 newform ----
mf7 = mfinit([7,3,-7],0); f7 = mfeigenbasis(mf7)[1]; L7 = lfunmf(mf7,f7);
Lf7_1 = lfun(L7,1); Lf7_2 = lfun(L7,2);
\\ ---- level 9 weight 3 newform: eta(3z)^6 ----
f9 = mffrometaquo(Mat([3,6]));
print("f9 params: ", mfparams(f9), "  coefs: ", mfcoefs(f9,13));
mf9 = mfinit(f9,0); L9 = lfunmf(mf9,f9);
Lf9_1 = lfun(L9,1); Lf9_2 = lfun(L9,2);
print("L(f7,1)=",Lf7_1); print("L(f7,2)=",Lf7_2);
print("L(f9,1)=",Lf9_1); print("L(f9,2)=",Lf9_2);
print("check FE f7: L1/L2 - sqrt7/(2Pi) = ", Lf7_1/Lf7_2 - sqrt(7)/(2*Pi));
print("check FE f9: L1/L2 - 3/(2Pi)     = ", Lf9_1/Lf9_2 - 3/(2*Pi));

L7chi = LK(-7,2); L3chi = LK(-3,2);
print("L(chi_-7,2)=",L7chi); print("L(chi_-3,2)=",L3chi);

nm7 = ["L(chi-7,2)","zeta(2)","Pi^2/sqrt7","Pi^3/sqrt7","Pi*L(chi-7,2)","L(f7,2)","L(f7,1)","Pi*L(f7,2)","Pi*L(f7,1)","Pi^2*L(f7,2)","sqrt7*L(f7,2)","Pi^2/(sqrt7)*L(chi-7,2)"];
bs7 = [L7chi, zeta(2), Pi^2/sqrt(7), Pi^3/sqrt(7), Pi*L7chi, Lf7_2, Lf7_1, Pi*Lf7_2, Pi*Lf7_1, Pi^2*Lf7_2, sqrt(7)*Lf7_2, Pi^2/sqrt(7)*L7chi];
nm3 = ["L(chi-3,2)","zeta(2)","Pi^2/sqrt3","Pi^3/sqrt3","Pi*L(chi-3,2)","L(f9,2)","L(f9,1)","Pi*L(f9,2)","Pi*L(f9,1)","Pi^2*L(f9,2)","sqrt3*L(f9,2)","Pi^2/sqrt3*L(chi-3,2)"];
bs3 = [L3chi, zeta(2), Pi^2/sqrt(3), Pi^3/sqrt(3), Pi*L3chi, Lf9_2, Lf9_1, Pi*Lf9_2, Pi*Lf9_1, Pi^2*Lf9_2, sqrt(3)*Lf9_2, Pi^2/sqrt(3)*L3chi];

two(tag, val, bs, nm) = {
  print("--- 2-term ratios for ", tag, " = ", val);
  for(i=1,#bs,
    my(r = val/bs[i], q = bestappr(r, 10^9));
    if(type(q)=="t_INT" || type(q)=="t_FRAC",
      my(res = abs(r - q));
      if(res < 1e-60, print("   HIT  ", tag, " = (", q, ") * ", nm[i], "   residual ", res),
         if(res < 1e-25, print("   weak ", tag, " = (", q, ") * ", nm[i], "   residual ", res)));
    );
  );
}
lin2(tag, val, bs, nm) = {
  print("--- lindep 2-term [val, b_i] for ", tag);
  for(i=1,#bs,
    my(v = lindep([val, bs[i]], 90));
    if(v != 0 && vecmax(abs(v)) < 10^8,
      my(res = abs(v[1]*val + v[2]*bs[i]));
      print("   ", nm[i], ": ", v, "  residual ", res));
  );
}
lin3(tag, val, bs, nm) = {
  print("--- lindep 3-term [val, b_i, b_j] for ", tag);
  for(i=1,#bs, for(j=i+1,#bs,
    my(v = lindep([val, bs[i], bs[j]], 80));
    if(v != 0 && v[1] != 0 && vecmax(abs(v)) < 10^7,
      my(res = abs(v[1]*val + v[2]*bs[i] + v[3]*bs[j]));
      if(res < 1e-70, print("   ", nm[i], " , ", nm[j], ": ", v, "  residual ", res)));
  ));
}
