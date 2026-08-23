default(parisizemax, 4000000000);
c4 = 864*t^5 + 40*t^4 - 256*t^3 + 96*t^2 - 16*t + 1;
c6 = -(5832*t^8 + 34560*t^7 - 30016*t^6 + 12624*t^5 - 4380*t^4 + 1280*t^3 - 240*t^2 + 24*t - 1);
A = -27*c4; B = -54*c6;
x0 = 144*t^2 - 72*t + 9;
print("x0 = ", factor(x0));
y2 = x0^3 + A*x0 + B;
print("y0^2 = ", factor(y2));
print("is square? ", issquare(y2));
r = 0; issquare(y2, &r); print("y0 = ", r);
/* also check 9-torsion possibility: is there a point P with 3P = the 3-torsion? */
/* and 6/12: no 2-torsion since cubic irreducible */
print();
print("psi3 second factor irreducible? ", factor(x^3 + (144*t^2 - 72*t + 9)*x^2 + (-46656*t^5 + 18576*t^4 - 6912*t^3 + 2592*t^2 - 432*t + 27)*x + (1259712*t^8 + 746496*t^7 - 449280*t^6 - 25920*t^5 + 92016*t^4 - 34560*t^3 + 6480*t^2 - 648*t + 27)));
quit;
