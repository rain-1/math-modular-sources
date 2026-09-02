default(realprecision,100);
E=ellinit([0,176,0,9216,147456]);
print("omega = ", E.omega);
Y=384*sqrt(-231);
P=[-384,Y];
print("on curve? ", ellisoncurve(E,P));
u=ellpointtoz(E,P);
print("u = ", u);
print("2u = ", 2*u);
print("u/omega1 = ", u/E.omega[1]);
print("u/omega2 = ", u/E.omega[2]);
\\ also x=1 for the E_m/F_m hosts is a branch point (2-torsion) - check
quit;
