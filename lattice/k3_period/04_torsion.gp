/* 04_torsion.gp -- Mordell-Weil torsion of the generic fibre over Q(t). */
default(parisizemax, 4000000000);
c4 = 864*t^5 + 40*t^4 - 256*t^3 + 96*t^2 - 16*t + 1;
c6 = -(5832*t^8 + 34560*t^7 - 30016*t^6 + 12624*t^5 - 4380*t^4 + 1280*t^3 - 240*t^2 + 24*t - 1);
dl = -t^4*(t-1)^6*(27*t^2-10*t+1)^3;
print("c4^3-c6^2-1728 Delta = ", c4^3-c6^2-1728*dl);
print("gcd(c4, 27t^2-10t+1) = ", gcd(c4, 27*t^2-10*t+1));
print("gcd(c4, t) = ", gcd(c4,t), "   gcd(c4,t-1) = ", gcd(c4,t-1));
A = -27*c4; B = -54*c6;
f = x^3 + A*x + B;
print();
print("factor of x^3+Ax+B over Q(t): ", factor(f));
ps3 = 3*x^4 + 6*A*x^2 + 12*B*x - A^2;
print();
print("psi_3 factors: ", factor(ps3));
/* 3-division: also try the 'x-coordinate of a 3-torsion' via the quartic roots */
print();
print("disc of psi3 wrt x: ", factor(polresultant(ps3, deriv(ps3,x), x)));
quit;
