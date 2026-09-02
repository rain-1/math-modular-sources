/* 12_cusptable.gp -- constant terms of Phi_D, Phi_new, Phi_new' at all four cusps. */
default(realprecision, 80);
G5 = znstar(5,1); ps4 = [G5,2]; ps4b = [G5,3];
P1 = mfeisenstein(3, ps4, 1); P2 = mfeisenstein(3, ps4b, 1);   /* my E_3^{1,psi4},E_3^{1,psi4bar} */
P3 = mfeisenstein(3, 1, ps4); P4 = mfeisenstein(3, 1, ps4b);   /* my E_3^{psi4,1},E_3^{psi4bar,1} */
mf1 = mfinit([5,3,ps4], 4); mf2 = mfinit([5,3,ps4b], 4);
cus = [[1,0;0,1], [2,1;5,3], [0,-1;1,0], [1,0;2,1]];
nam = ["oo ", "2/5", "0  ", "1/2"];
s5 = sqrt(5); ph5 = (11+5*s5)/2; phm5 = (5*s5-11)/2;
c0(mf, F, g) = my(p=0); my(e = mfslashexpansion(mf, F, g, 1, 0, &p)); [subst(lift(e[1]),t,I)*1.0, p[2]];
{ for(j=1,4, my(g=cus[j]);
    my(a1=c0(mf1,P1,g), a2=c0(mf2,P2,g), a3=c0(mf1,P3,g), a4=c0(mf2,P4,g));
    my(w=a1[2]);
    my(cD  = (1/2+I)*a1[1] + (1/2-I)*a2[1]);
    my(cN  = (1+I*ph5)*a3[1] + (1-I*ph5)*a4[1]);
    my(cNp = (1-I*phm5)*a3[1] + (1+I*phm5)*a4[1]);
    print("cusp ", nam[j], " (width ", w, "):");
    print("    c0(E_3^{1,psi4})    = ", a1[1]);
    print("    c0(E_3^{1,psi4bar}) = ", a2[1]);
    print("    c0(E_3^{psi4,1})    = ", a3[1]);
    print("    c0(E_3^{psi4bar,1}) = ", a4[1]);
    print("    c0(Phi_D)     = ", cD);
    print("    c0(Phi_new)   = ", cN);
    print("    c0(Phi_new')  = ", cNp);
  ); }
print("");
print("rho = 1/(125*5^(3/4)*phi^(5/2)) = ", 1/(125*5^(3/4)*ph5^(1/2)));
print("check c0(P3,cusp0) = rho*(-1+I*phi^5) = ", (1/(125*5^(3/4)*ph5^(1/2)))*(-1+I*ph5));
quit;
