/* b13_conifold.gp -- xi_infty(207) against the Gamma-class / conifold basis at the
   dominant conifold z_- = (349-85 sqrt17)/2^17 (exponents 0,1,1,2), over Q and Q(sqrt17). */
default(realprecision, 400);
xi0 = eval(Str(read("b2_xi_arch.txt")));
default(realprecision, 335);
xi = xi0*1.0;
s17 = sqrt(17.);
zm = (349-85*s17)/131072;  zp = (349+85*s17)/131072;  zc = 1/53248.;
Lm = log(abs(1/zm)); Lp = log(1/zp); Lc = log(1/zc);
print("Lm=",Lm,"  Lp=",Lp,"  Lc=",Lc);
z2 = zeta(2); z3 = zeta(3); P=Pi;
try(nm, v) = { my(r = lindep(v, 300));
   if(type(r)=="t_COL" && r[1]!=0 && vecmax(abs(r)) < 10^(floor(300/(#v-1))-2),
      print("  ", nm, " : ", r~), print("  ", nm, " : -")); }

print("\n=== Gamma-class basis at z_- ===");
try("[1,Lm,Lm^2,Lm^3,z2,z2*Lm,z3]", [xi,1,Lm,Lm^2,Lm^3,z2,z2*Lm,z3]);
try("[1,Lm,Lm^2,Lm^3,P^2,P^2*Lm,z3]", [xi,1,Lm,Lm^2,Lm^3,P^2,P^2*Lm,z3]);
try("[1,Lm,Lm^2,Lm^3,z2*Lm,z3]",      [xi,1,Lm,Lm^2,Lm^3,z2*Lm,z3]);
try("[1,Lm,Lm^2,Lm^3,z3]",            [xi,1,Lm,Lm^2,Lm^3,z3]);
try("[1,Lm,Lm^2,Lm^3]",               [xi,1,Lm,Lm^2,Lm^3]);
try("[1,z2,z3]",                      [xi,1,z2,z3]);
try("[1,Lm,z2,z3]",                   [xi,1,Lm,z2,z3]);
print("\n=== same, with sqrt17 twists ===");
try("[1,s,Lm,s*Lm,Lm^2,s*Lm^2,z3,s*z3]", [xi,1,s17,Lm,s17*Lm,Lm^2,s17*Lm^2,z3,s17*z3]);
try("[1,s,z2,s*z2,z3,s*z3]",             [xi,1,s17,z2,s17*z2,z3,s17*z3]);
try("[1,s,Lm,s*Lm,z3,s*z3]",             [xi,1,s17,Lm,s17*Lm,z3,s17*z3]);
try("[1,s,Lm,s*Lm,Lm^2,s*Lm^2,Lm^3,s*Lm^3]", [xi,1,s17,Lm,s17*Lm,Lm^2,s17*Lm^2,Lm^3,s17*Lm^3]);
print("\n=== using log of the algebraic singular value beta=85sqrt17-349 ===");
{ my(bb = log(85*s17-349), aa = log(349+85*s17), l2=log(2.), l13=log(13.));
  print("  log|beta|=",bb,"  log alpha=",aa);
  try("[1,bb,bb^2,bb^3,z3]", [xi,1,bb,bb^2,bb^3,z3]);
  try("[1,l2,bb,z3,z2]",     [xi,1,l2,bb,z3,z2]);
  try("[1,l2,l13,bb,z3,z2]", [xi,1,l2,l13,bb,z3,z2]);
  try("[1,l2,l13,aa,z3,z2,P^2]", [xi,1,l2,l13,aa,z3,z2,P^2]);
}
print("\n=== xi times powers / with Pi^3 ===");
try("[xi*P^3 ; 1,Lm,Lm^2,Lm^3,z3]", [xi*P^3,1,Lm,Lm^2,Lm^3,z3]);
try("[xi/z3 ; 1,Lm,Lm^2,Lm^3]",     [xi/z3,1,Lm,Lm^2,Lm^3]);
try("[xi ; z3,P^2*Lm,Lm^3,1]",      [xi,z3,P^2*Lm,Lm^3,1]);
print("\n=== chi/H^3 = -120/13 (AESZ 99) tests ===");
print("  xi/(zeta(3)) = ", xi/z3, "   *13/120 = ", xi/z3*13/120, "   *(-13/120)=",-xi/z3*13/120);
print("  bestappr(xi/zeta(3)*13,10^12) = ", bestappr(xi/z3*13,10^12));
quit
