/* 10_cusps.gp -- constant terms at the four cusps of Gamma_1(5), weight 3.
   NOTE PARI's mfeisenstein(k,C1,C2) has c_n = sum_{d|n} C1(d)C2(n/d)d^{k-1},
   i.e. the OPPOSITE order to the eis(ch1,ch2) used in the hostscan scripts.
   my P1 = E_3^{1,psi4}    = mfeisenstein(3, ps4,  1)
   my P2 = E_3^{1,psi4bar} = mfeisenstein(3, ps4b, 1)
   my P3 = E_3^{psi4,1}    = mfeisenstein(3, 1, ps4 )
   my P4 = E_3^{psi4bar,1} = mfeisenstein(3, 1, ps4b)                        */
default(realprecision, 38);
G5 = znstar(5,1);
ps4 = [G5, 2]; ps4b = [G5, 3];
P1 = mfeisenstein(3, ps4, 1); P2 = mfeisenstein(3, ps4b, 1);
P3 = mfeisenstein(3, 1, ps4); P4 = mfeisenstein(3, 1, ps4b);
print("P1 = E_3^{1,psi4}  q-exp : ", mfcoefs(P1, 6));
print("P3 = E_3^{psi4,1}  q-exp : ", mfcoefs(P3, 6));
mf1 = mfinit([5,3,ps4], 4); mf2 = mfinit([5,3,ps4b], 4);
cus = [[1,0;0,1], [2,1;5,3], [0,-1;1,0], [1,0;2,1]];
nam = ["oo=1/0", "2/5   ", "0=0/1 ", "1/2   "];
c0f(mf, F, g) = my(p=0); my(e = mfslashexpansion(mf, F, g, 3, 0, &p)); [e, p];
{ for(j=1, 4,
   my(g = cus[j]);
   my(r1 = c0f(mf1,P1,g), r2 = c0f(mf2,P2,g), r3 = c0f(mf1,P3,g), r4 = c0f(mf2,P4,g));
   print("cusp ", nam[j], "   params(P1,P2,P3,P4) = ", r1[2], r2[2], r3[2], r4[2]);
   print("   P1|g : ", r1[1]);
   print("   P2|g : ", r2[1]);
   print("   P3|g : ", r3[1]);
   print("   P4|g : ", r4[1]);
 ); }
quit;
