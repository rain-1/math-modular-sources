print("=== C6: E: v^2=u^3-1, z=u^3 ===");
e = ellinit([0,0,0,0,-1]);
print("j(E) = ", e.j, "   disc = ", e.disc, "   c4 = ", e.c4, " c6 = ", e.c6);
print("-- ramification of z=u^3 on E (degree 6) --");
print("over z=0: u=0 => v^2=-1 => v=+-I, two points; u is a uniformiser there (v!=0), z=u^3 => e=3 each.");
print("   check: u^3-1 at u=0 is ", subst('u^3-1,'u,0), " (nonzero) so v!=0. #pts = ", 2);
print("over z=1: u^3=1 => v^2=0 => v=0, Weierstrass pts; roots of u^3-1: ", polroots('u^3-1));
print("   local: z-1 = u^3-1 = v^2, so e=2 at each of 3 points.");
print("over z=oo: single point at infinity, ord(u)=-2, ord(v)=-3, ord(z)=ord(u^3)=-6 => e=6.");
{
my(rh = 6*(2*0-2) + (2*(3-1) + 3*(2-1) + 1*(6-1)));
print("Riemann-Hurwitz: 2g-2 = 6*(-2) + [2*2 + 3*1 + 1*5] = ", rh, "  => g = ", (rh+2)/2);
}
print("sum of degrees: over 0: 3+3=6 ; over 1: 2+2+2=6 ; over oo: 6 -> all = deg 6  OK");

print();
print("=== C8: Belyi phi(x) = -x^6/(2x^3+1) ===");
num = -'x^6; den = 2*'x^3+1;
print("deg phi = max(", poldegree(num), ",", poldegree(den), ") = 6");
print("phi zeros: factor(num) = ", factor(num));
print("phi poles (finite): factor(den) = ", factor(den), "  ; at x=oo ord(phi) = ", poldegree(den)-poldegree(num), " (pole of order 3)");
print("phi-1 numerator = ", factor(num-den), "   [= -(x^3+1)^2]");
print("check -(x^3+1)^2 == num-den : ", num-den == -('x^3+1)^2);
print("factor(x^3+1) = ", factor('x^3+1));
print("=> fibre over 0 : x=0 with mult 6                -> cycle type (6)");
print("=> fibre over 1 : 3 roots of x^3+1, each mult 2   -> cycle type (2,2,2)");
print("=> fibre over oo: 3 roots of 2x^3+1 (simple) + x=oo with mult 3 -> (3,1,1,1)");
{ my(rh = 6*(-2) + ((6-1) + 3*(2-1) + (3-1)));
print("RH for phi: 2g-2 = ", rh, " => g = ", (rh+2)/2, " (genus 0, consistent)"); }

print();
print("=== C10: Lyons factorisation ===");
aa = 1-6*'T^2-3*'T^4; bb = 16*'T^6;
lhs = aa^2+12*bb;
rhs = (3*'T^2+1)*(3*'T^6+75*'T^4-15*'T^2+1);
print("A^2+12B = ", lhs);
print("RHS      = ", rhs);
print("EQUAL? ", lhs==rhs, "   difference = ", lhs-rhs);
print("factor(A^2+12B) over Q[T] = ", factor(lhs));
print("-- 2-isogenous curve hatE: y^2 = x^3 -2A x^2 + (A^2-4B) x --");
{
my(a2=-2*aa, a4=aa^2-4*bb, b2, b4, c4, disc);
b2 = 4*a2; b4 = 2*a4;
c4 = b2^2-24*b4;
print("b2 = ", b2, "  b4 = ", b4);
print("c4(hatE) = ", c4);
print("c4 / (A^2+12B) = ", c4/lhs, "   -> proportional? ", type(c4/lhs)=="t_INT" || poldegree(c4/lhs)==0);
disc = 16*a4^2*(a2^2-4*a4);
print("Disc(hatE) = ", factor(disc));
print("c4(E_tau)  = ", (4*aa)^2-24*(2*bb), " = 16*(A^2-3B)");
}
print("-- specialise T^2 = -1/3 --");
{
my(s=-1/3, A0, B0, c4h, dh, c4e, de);
A0 = 1-6*s-3*s^2; B0 = 16*s^3;
print("A = ", A0, "  B = ", B0);
c4h = 16*(A0^2+12*B0);
dh  = 16*(A0^2-4*B0)^2*(4*A0^2-4*(A0^2-4*B0));
print("c4(hatE) = ", c4h, "   Disc(hatE) = ", dh, "  nonzero? ", dh!=0);
print("j(hatE) = c4^3/Disc = ", c4h^3/dh);
c4e = 16*(A0^2-3*B0);
de  = 16*B0^2*(A0^2-4*B0);
print("c4(E_tau) = ", c4e, "  Disc(E_tau) = ", de, "  j(E_tau) = ", c4e^3/de, " = ", c4e^3/de*1.0);
}
quit();
