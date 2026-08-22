default(parisizemax, 4000000000);
default(realprecision, 60);
\\ minimal polynomial over Q of lambda (root of lam^2-(248+110 sqrt5)lam+1)
mp = (x^2-248*x+1)^2 - 60500*x^2;
print("min poly over Q of lambda : ", mp);
print("  irreducible over Q ? ", polisirreducible(mp));
print("  roots : ", polroots(mp));
print("  constant term (= N_{K(lam)/Q}(lambda)) = ", polcoeff(mp,0), "   -> lambda is a UNIT of absolute norm 1");
print("  poldisc = ", factor(poldisc(mp)));
K = nfinit(y^2-y-1);
disc = Mod(122000+54560*(2*y-1), y^2-y-1);
print("char-poly discriminant (248+110sqrt5)^2-4 = 122000+54560 sqrt5 ; in basis [1,w] : ", lift(disc));
print("  is it a square in K ?  nfroots(K, x^2 - disc) = ", nfroots(K, x^2 - (122000+54560*(2*y-1))));
nf4 = nfinit(mp);
print("quartic field K(lambda): degree ", poldegree(mp), ", signature ", nf4.sign, ", disc = ", nf4.disc, " = ", factor(nf4.disc));
print("  contains sqrt5 ? ", #nfroots(nf4, x^2-5) > 0);
print("  Galois group of mp : ", polgalois(mp));
\\ place-by-place products
default(realprecision,60);
r = polroots(mp);
print("all four conjugates of lambda : ", vecsort(vector(4,i,real(r[i]))));
l1v1 = 493.965453091904187843965037193676477066833035760879891619376;
l2v1 = 1/l1v1;
l1v2 = 1.19733277825870061570520501991715891206360693969050568280042;
l2v2 = 1/l1v2;
print("small root at v1 = ", l2v1, "  small root at v2 = ", l2v2);
print("place-by-place product of the SMALL roots = ", l2v1*l2v2, "  its inverse = ", 1/(l2v1*l2v2));
print("  |N|^(1/2) (place-by-place, small roots) = ", sqrt(l2v1*l2v2), "  inverse = ", 1/sqrt(l2v1*l2v2));
print("full absolute norm over the quartic field = ", l1v1*l2v1*l1v2*l2v2);
