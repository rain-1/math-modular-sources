/* Exact location of the nearest archimedean singularity for each curve:
   elliptic points (j=0 or 1728 with ramification index 1) resp. non-cuspidal
   values of the hauptmodul at the other cusps.  x = (Delta(N tau)/Delta(tau))^{1/(N-1)}.
   j-formulas in terms of x = 1/t, t the usual eta quotient hauptmodul.        */
default(realprecision, 40);
print("X_0(2):  j = (256x+1)^3/x ;  elliptic pts of order 2 (j=1728):");
f2 = (256*x+1)^3 - 1728*x;
print("   factor: ", factor(f2));
print("   roots : ", polroots(f2), "   |.| = ", [abs(z) | z <- polroots(f2)], "   2^-6 = ", 1.0/64);
print("");
print("X_0(3):  j = (27x+1)(243x+1)^3/x ; j=1728:");
f3 = (27*x+1)*(243*x+1)^3 - 1728*x;
print("   factor: ", factor(f3));
print("   roots : ", polroots(f3), "   |.| = ", [abs(z) | z <- polroots(f3)], "   3^-3 = ", 1.0/27);
print("");
print("X_0(5):  j = (3125x^2+250x+1)^3/x ; j=1728:");
f5 = (3125*x^2+250*x+1)^3 - 1728*x;
print("   factor: ", factor(f5));
print("   |roots| = ", vecsort([abs(z) | z <- polroots(f5)]), "   5^{-3/2} = ", 5.0^(-1.5));
print("");
print("X_0(7):  j = (49x^2+13x+1)(2401x^2+245x+1)^3/x ; j=0:");
f7a = 49*x^2+13*x+1; f7b = 2401*x^2+245*x+1;
print("   elliptic (e=1) factor 49x^2+13x+1 roots: ", polroots(f7a), " |.| = ", [abs(z) | z <- polroots(f7a)], "  7^-1 = ", 1.0/7);
print("   non-elliptic (e=3) factor roots: ", polroots(f7b), " |.| = ", [abs(z) | z <- polroots(f7b)]);
quit;
