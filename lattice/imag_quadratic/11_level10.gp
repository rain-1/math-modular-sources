\\ 11_level10.gp -- the first Q(i)-rational fold-regular source with a NON-REAL period.
\\ Level 10, weight 3, nebentypus psi_4 (mod 5, odd, order 4).  Inner orientation with
\\ the oldform shift:   Phi = E_3^{psi_4,1}(tau) - 8 E_3^{psi_4,1}(2 tau),
\\ coefficients in Z[i]; a_0 = 0 at infinity (psi_1 != 1) and at the cusp 0
\\ (sum_d c_d d^{-3} = 1 - 8/8 = 0).  Period at the near cusp: -1/2 L(2,psi_4), non-real.
\\ Exactly this object at level 5 does not exist (04_g15.out: no oldform shift available),
\\ and at level 10 there is no four-point genus-zero host (01_groups.out).
default(realprecision,60);
G5 = znstar(5,1);
ps4 = znconreylog(G5,2);
G1 = znstar(1,1);
tr = znconreylog(G1,1);
Ein = mfeisenstein(3,[G1,tr],[G5,ps4]);
Eout = mfeisenstein(3,[G5,ps4],[G1,tr]);
print("INNER  E_3^{psi_4,1}, c(m) = sum_{d|m} psi_4(m/d) d^2 : ", mfcoefs(Ein,8));
print("OUTER  E_3^{1,psi_4}, c(m) = sum_{d|m} psi_4(d) d^2   : ", mfcoefs(Eout,8));
G10 = znstar(10,1);
mf10 = mfinit([10,3,[G10,znconreylog(G10,7)]],4);
print("dim M_3^Eis(Gamma_1(10), psi_4) = ", mfdim(mf10));
Phi = mflinear([Ein,mfbd(Ein,2)],[1,-8]);
print();
print("Phi = E_in(tau) - 8 E_in(2 tau) : ", mfcoefs(Phi,10));
print("   all coefficients in Z[i]      : yes (visible above)");
print("   in the level-10 psi_4 space   : ", mftobasis(mf10,Phi));
print("   a_0 at infinity               : ", mfcoefs(Phi,0));
pp = 0;
v = mfslashexpansion(mf10, Phi, [0,-1;1,0], 2, 1, &pp);
print("   a_0 at the cusp 0             : ", v[1]);
print("   => FOLD-REGULAR, and Q(i)-rational.");
print();
Psi = mflinear([Eout,mfbd(Eout,2)],[1,-1]);
print("(companion outer direction E_out(tau) - E_out(2 tau): a_0 at oo = ",mfcoefs(Psi,0));
pp = 0;
v2 = mfslashexpansion(mf10, Psi, [0,-1;1,0], 2, 1, &pp);
print(" a_0 at cusp 0 = ", v2[1], "  -- also fold-regular, but its period is zeta(2) L(0,psi_4), i.e. zeta(2) times an algebraic number)");
print();
L2 = lfun(lfuncreate([G5,ps4]),2);
print("L(2,psi_4) = ", L2);
print("Re = ", real(L2));
print("Im = ", imag(L2), "   -- NOT real.");
print("A K-rational Apery system carrying it would prove L(2,psi_4) not in Q(i),");
print("a statement that no real-period system can give (K cap R = Q).  The four-point");
print("geometry needed to run it stops at N = 9.");
quit;
