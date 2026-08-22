/* value of the X_0(9) hauptmodul u = (eta(9 tau)/eta(tau))^3 at the cusps 1/3, 2/3,
   and at the Fricke fixed point tau = i/3.   Also X_0(4): z at cusp 1/2. */
default(realprecision, 40);
uu(t) = (eta(9*t, 1)/eta(t, 1))^3;
zz(t) = (eta(4*t, 1)/eta(t, 1))^8;
print("X_0(9):  u(i/3) [Fricke fixed pt] = ", uu(I/3), "   3^{-3/2} = ", 3.0^(-1.5));
for(e = 1, 5, my(t = 10.0^(-e)); print("  u(1/3 + ", t, "i) = ", uu(1/3 + t*I)));
for(e = 1, 5, my(t = 10.0^(-e)); print("  u(2/3 + ", t, "i) = ", uu(2/3 + t*I)));
print("  |limit| vs 3^{-3/2} = ", 3.0^(-1.5), " ;  arg/pi of limit ~ +-2/3 ?");
print("");
print("X_0(4):  z(i/2) [Fricke fixed pt] = ", zz(I/2), "   2^{-4} = ", 1.0/16);
for(e = 1, 5, my(t = 10.0^(-e)); print("  z(1/2 + ", t, "i) = ", zz(1/2 + t*I)));
quit;
