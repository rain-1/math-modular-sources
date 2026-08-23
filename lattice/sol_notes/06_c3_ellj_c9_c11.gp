print("=== C6 (j-invariant) ===");
e = ellinit([0,0,0,0,-1]);
print("E: y^2 = x^3 - 1 :  c4 = ", e.c4, "  c6 = ", e.c6, "  disc = ", e.disc, "  j = ", e.j);
e2 = ellinit([0,0,0,0,16]);  \\ sanity, another j=0 curve
print("sanity j(y^2=x^3+16) = ", e2.j);

print();
print("=== C9: v^2 = u^3 - 1 in Q(x,y)/(y^6+2x^3+1) ===");
\\ work in Q[x,y]/(y^6 + 2x^3 + 1); clear denominators: check (x^3+1)^2 * y^6 ?= (x^6 - y^6) * y^6 ... do it properly
\\ v^2 - (u^3-1) = (x^3+1)^2/y^6 - x^6/y^6 + 1 = ((x^3+1)^2 - x^6 + y^6)/y^6
numexpr = ('x^3+1)^2 - 'x^6 + 'y^6;
print("numerator N = (x^3+1)^2 - x^6 + y^6 = ", numexpr);
rel = 'y^6 + 2*'x^3 + 1;
print("relation R = y^6+2x^3+1 = ", rel);
print("N - R = ", numexpr - rel, "   => N = R exactly? ", numexpr == rel);
print("Hence N == 0 in the function field, so v^2 = u^3 - 1 EXACTLY.");
\\ independent: reduce N modulo R as polynomial in y
print("N mod R (as poly in y) = ", numexpr % rel);
\\ also check z = u^3 = phi(x)
print("u^3 = x^6/y^6 and y^6 = -(2x^3+1)  =>  u^3 = -x^6/(2x^3+1) = phi(x).");

print();
print("=== C11: numerical invariants ===");
{
my(lam=1, d0=10, d1=0, gB=1, gF=2, chi, ee, k2, b1, b2, h11, pg, q);
chi = lam + (gB-1)*(gF-1);
print("chi(O_X) = lambda + (g(B)-1)(g(F)-1) = ", lam, " + ", (gB-1)*(gF-1), " = ", chi, "   [claim 1] ", chi==1);
ee = 4*(gB-1)*(gF-1) + d0 + d1;
print("e(X) = 4(g_B-1)(g_F-1) + delta = ", ee, "   [claim 10] ", ee==10);
k2 = 12*chi - ee;
print("K^2 = 12chi - c2 = ", k2, "   [claim 2] ", k2==2);
q = gB; pg = chi + q - 1;
print("q = ", q, "  p_g = chi+q-1 = ", pg, "   [claim 1,1] ", q==1 && pg==1);
b1 = 2*q; print("b1 = 2q = ", b1, "   [claim 2] ", b1==2);
b2 = ee - 2 + 2*b1;
print("b2 = e - 2 + 2 b1 = ", b2, "   [claim 12] ", b2==12);
h11 = b2 - 2*pg;
print("h^{1,1} = b2 - 2 p_g = ", h11, "   [claim 10] ", h11==10);
print("Mumford: 10*lambda = delta0 + 2*delta1 : ", 10*lam, " =? ", d0+2*d1, "  -> ", 10*lam == d0+2*d1);
print("Noether: chi = (K^2+c2)/12 = ", (k2+ee)/12);
print("h^{1,1} vs rho<=h^{1,1}: b2 = 2 p_g + h11 = ", 2*pg+h11);
}
print();
print("=== C12: dual graph Betti numbers vs rank N ===");
{
print("theta graph: V=2, E=3 -> b1 = E-V+1 = ", 3-2+1, "  vs rank(N0)=2 -> ", (3-2+1)==2);
print("4-cycle    : V=4, E=4 -> b1 = E-V+1 = ", 4-4+1, "  vs rank(Ninf)=1 -> ", (4-4+1)==1);
print("total nodes = 3+3+4 = ", 3+3+4, " = e(X) = 10 ? ", 10==10);
print("component excesses: (2-1)+(2-1)+(4-1) = ", 1+1+3);
print("Q(-1)^(2+5)=Q(-1)^7, b2=12 -> dim H_par = ", 12-7);
}
quit();
