default(realprecision, 200);
XI = 0.679442502761532499455298305916444856664282472334996092997440915969317355781462466592540405781946691328149574580693860939;
G = Catalan; Z2 = zeta(2); L2 = log(2);
{try(nm, v, h) = my(r = lindep(v, 110));
  print("  ", nm, " -> ", r~, "  |sum| = ", abs(sum(i=1,#v, r[i]*v[i])));}
try("[1,xi,G,z2]", [1., XI, G, Z2]);
try("[1,xi,G,Pi*L2]", [1., XI, G, Pi*L2]);
try("[1,xi,G,Pi*L2,L2^2]", [1., XI, G, Pi*L2, L2^2]);
try("[1,xi,G,Pi*L2,Pi^2,L2^2,Pi,L2]", [1., XI, G, Pi*L2, Pi^2, L2^2, Pi, L2]);
try("[xi,G,Pi*L2]", [XI, G, Pi*L2]);
try("[xi,G,Pi]", [XI, G, Pi]);
try("[xi,G,L2]", [XI, G, L2]);
try("[xi,Pi,L2,1]", [XI, Pi, L2, 1.]);
try("[xi,G/Pi,1,Pi]",[XI, G/Pi, 1., Pi]);
print("xi*Pi = ", XI*Pi); print("xi/Pi = ", XI/Pi);
print("xi*4 = ", XI*4); print("xi-log2 = ", XI-L2);
print("xi vs G/2+something: xi-G/2 = ", XI-G/2, "   xi+G/2 = ", XI+G/2);
try("[1,xi,G,z2,Pi*L2,L2^2,zeta(3)/Pi]", [1., XI, G, Z2, Pi*L2, L2^2, zeta(3)/Pi]);
quit;
