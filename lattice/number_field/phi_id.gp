default(parisizemax,4000000000);
default(realprecision,60);
w=Mod(y,y^2-y-1); ph=w;
print("4*phi^10 + 2 = ", lift(4*ph^10+2), "   and 248+110sqrt5 in basis[1,w] = ", lift(248+110*(2*w-1)));
print("EQUAL ? ", 4*ph^10+2 == 248+110*(2*w-1));
p5=((1+sqrt(5))/2)^5;
print("(phi^5 + sqrt(phi^10+1))^2 = ", (p5+sqrt(p5^2+1))^2);
print("lambda_1 at v1            = 493.965453091904187843965037193676477066833035760879891619376");
print("sqrt(lambda_1) - 1/sqrt(lambda_1) = ", sqrt(493.965453091904187843965037193676477066833035760879891619376)-1/sqrt(493.965453091904187843965037193676477066833035760879891619376), "   2*phi^5 = ", 2*p5);
