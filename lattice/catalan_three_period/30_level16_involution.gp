/* 30_level16_involution.gp   (2026-08-23)
   Independent check of the level-16 descending involution sigma(x) = -x/(4x+1)
   and of the geometry of the symmetrised y-line y = 4x^2/(4x+1).

   [1] sigma is induced by  tau -> tau + 1/2.
   [2] cusp values of the two hosts, computed at every cusp (eta(z,1)):
         level 8  t8  : {0, 1/8, 1/4, oo}      -> 2s = 1/2 is NOT a cusp (free)
         level 16 x16 : {0,-1/4,-1/2,(-1+-i)/4,oo} -> 2s = -1/2 IS a cusp
   [3] y = 4x^2/(4x+1) is a series in q^2 and y(tau) = -4*t8(2tau+1/2):
         the level-16 y-line is the level-8 x-line rescaled by -4, i.e. Y_0(8):
         4 punctures {0,-1/2,-1,oo}, NO orbifold point.
   [4] level-8 sigma_8(x) = x/(4x-1) is not modular (sigma_8(1/8) = -1/4).
*/
\p 40
NT = 60;  q = 'q;
et(k) = eta(q^k + O(q^(NT+2)));
x   = q*et(2)*et(16)^2/(et(1)^2*et(8)) + O(q^NT);
t8  = q*et(1)^4*et(4)^2*et(8)^4/et(2)^10 + O(q^NT);
print("x   = ", x + O(q^12));
print("t8  = ", t8 + O(q^12));

print("\n[1] sigma = T^(1/2) ?");
xhalf = subst(x, q, -q);
print("    x(tau+1/2) + x/(4x+1) = ", xhalf + x/(4*x+1));
print("    (zero to O(q^", NT, ")  =>  sigma(x) = -x/(4x+1) is tau -> tau+1/2)");
print("    fixed points of sigma : roots of 4x^2+2x : ", polroots(4*'x^2+2*'x)~);
print("    fixed points of sigma_8(x)=x/(4x-1) : ", polroots(4*'x^2-2*'x)~);

print("\n[2] cusp values");
E(z) = eta(z,1);
T8(z)  = E(z)^4*E(4*z)^2*E(8*z)^4/E(2*z)^10;
X16(z) = E(2*z)*E(16*z)^2/(E(z)^2*E(8*z));
cv(f,a,c,T) = my(g=bezout(a,c), d=g[1], b=-g[2]); f((a*I*T+b)/(c*I*T+d));
print("    level 8, Gamma_0(8), cusps 0,1/2,1/4,oo :");
foreach([[0,1],[1,2],[1,4],[1,8]], P, print("       ",P[1],"/",P[2]," -> ",cv(T8,P[1],P[2],24)));
print("    level 16, Gamma_0(16), cusps 0,1/2,1/4,3/4,1/8,oo :");
foreach([[0,1],[1,2],[1,4],[3,4],[1,8],[1,16]], P, print("       ",P[1],"/",P[2]," -> ",cv(X16,P[1],P[2],24)));

print("\n[3] the symmetrised coordinate");
y = 4*x^2/(4*x+1) + O(q^NT);
print("    y = ", y + O(q^14));
print("    odd part of y (must be 0): ", sum(n=1,NT-1, if(n%2==1, polcoeff(y,n)*q^n, 0)));
u = subst(subst(t8,q,-q), q, q^2) + O(q^NT);   /* t8(2tau+1/2) */
print("    y + 4*t8(2tau+1/2) = ", y + 4*u);
print("    => y(tau) = -4*t8(2tau+1/2); y = 4q^2 + ..., canonical cusp");
print("       parameter Q = q^2 (width 1/2), so |dy/dQ|(0) = 4 exactly.");

print("\n[4] level-8 involution is not modular");
print("    sigma_8(1/8) = ", (1/8)/(4/8-1), "  not a cusp value of t8");
print("    t8(tau+1/2) - t8/(4*t8-1) = ", subst(t8,q,-q) - t8/(4*t8-1) + O(q^8));
quit
