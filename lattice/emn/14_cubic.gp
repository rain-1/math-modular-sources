/* (j) sec.13: Ktil(z) = H(R3(z)) + 3H(z),  R3(z) = z(3-2z)^2.
   Names: hser, r3, ktil (avoid H, K, S). */
default(parisizemax, 12*10^9);
NORD = 420;
amom(n) = sum(k = 0, n, binomial(n,k)/(2*k+1)) / (2^(n+1)*(n+1));
av = vector(NORD+2); for(n = 0, NORD+1, av[n+1] = amom(n));

hser = sum(n = 0, NORD, av[n+1]*'z^(n+1)) + O('z^(NORD+2));
r3 = 'z*(3-2*'z)^2;
print("R3(z) = ", r3, ";  R3(1) = ", subst(r3,'z,1), ";  R3'(1) = ", subst(deriv(r3),'z,1));
print("R3(z)-1 factors: ", factor(4*'z^3-12*'z^2+9*'z-1));
print("R3(z)-2 factors: ", factor(4*'z^3-12*'z^2+9*'z-2));
gettime();
comp = subst(hser, 'z, r3 + O('z^(NORD+2)));
ktil = comp + 3*hser;
print("composition time: ", gettime(), " ms");

print();
print("=== (j1) coefficients rational for n>=1; v_2 slope ===");
{
my(vv = vector(NORD));
for(n = 1, NORD, vv[n] = valuation(polcoef(ktil, n, 'z), 2));
print("  n, v_2([z^n] Ktil), v_2 - n:");
for(i = 1, 20, print("   n=", i, "  v_2=", vv[i], "  v_2-n=", vv[i]-i));
foreach([25,50,75,100,150,200,250,300,350,400], n, print("   n=", n, "  v_2=", vv[n], "  v_2-n=", vv[n]-n, "  v_2/n=", 1.0*vv[n]/n));
my(mn = 10^9, mx = -10^9);
for(n = 1, NORD, my(d = vv[n]-n); if(d<mn, mn=d); if(d>mx, mx=d));
print("  range of v_2 - n over 1<=n<=", NORD, ": [", mn, ", ", mx, "]");
/* is the deviation O(log n)? */
print("  v_2(n) contributions: compare v_2 - n with -v_2(n!) etc.");
for(i = 1, 8, print("   n=", 2^i, "  v_2-n = ", vv[2^i]-2^i));
KV = vv;
}

print();
print("=== (j2) compare with H itself and with the pure part ===");
{
for(i = 1, 10, print("   n=", i, "  v_2([z^n]H) = ", valuation(polcoef(hser,i,'z),2), "   v_2([z^n]H(R3)) = ", valuation(polcoef(comp,i,'z),2)));
}

print();
print("=== (j3) radius of convergence / nearest singularity ===");
{
print("  R3^{-1}(1) = roots of 4z^3-12z^2+9z-1: ", polroots(4*'z^3-12*'z^2+9*'z-1));
print("  1 - sqrt(3)/2 = ", 1 - sqrt(3)/2, ";  1 + sqrt(3)/2 = ", 1 + sqrt(3)/2);
print("  ratio test on Ktil coefficients (|c_n|^{-1/n}):");
foreach([50,100,200,300,400], n, my(c = polcoef(ktil,n,'z)); print("   n=", n, "  |c_n|^{-1/n} = ", exp(-log(abs(c*1.0))/n)));
}

print();
print("=== (j4) rigidity: is there a nonconstant combination killing all three branches? ===");
print("  the cubic fibre over R3(z) is {z, z1, z2}; the note claims the only");
print("  coefficient pattern annihilating all of 1,1-sqrt3/2,1+sqrt3/2 is (1,3,3,3)");
print("  i.e. the full trace, which equals the constant 10G.");
write("/home/ubuntu/code/math-modular-sources/lattice/emn/ktil_v2.txt", KV);
quit;
