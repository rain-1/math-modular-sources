/* the singular set of H: is z=1 the ONLY finite singularity? */
default(parisizemax, 4*10^9);
amom(n) = sum(k = 0, n, binomial(n,k)/(2*k+1))/(2^(n+1)*(n+1));
print("=== asymptotics of a_n: a single dominant singularity at z=1 (no oscillation) ===");
foreach([10,50,100,200,400,800], n, print("   n=", n, "   n^2 a_n = ", 1.0*n^2*amom(n), "   (2 n^2 a_n = ", 2.0*n^2*amom(n), ")"));
print("   -> n^2 a_n -> 1/2 monotonically: the only singularity on |z|=1 is z=1,");
print("      and the local type is (1-z)log(1-z) (coefficients ~ 1/(2n^2), no alternation).");
print();
print("=== homogeneous solutions of z(z-2)Y''+(z-1)Y'=0 ===");
{
/* Y = uni(z) = 2 asin(sqrt(z/2)) : substitute z = 2*w^2 so uni = 2 asin(w) is a series */
my(w = 'w + O('w^30));
my(uni = 2*asin(w));       /* = uni(z) with z = 2 w^2 */
/* d/dz = (1/(4w)) d/dw ;  Y'' = (1/(4w)) d/dw [ (1/(4w)) dY/dw ] */
my(d1 = deriv(uni,'w)/(4*'w));
my(d2 = deriv(d1,'w)/(4*'w));
my(zz = 2*'w^2);
my(res = zz*(zz-2)*d2 + (zz-1)*d1);
print("   residual of the homogeneous ODE on uni = 2 asin(w), z = 2w^2 : ", res);
}
print();
print("=== H is regular at z=2: H(z) = 2G - H(2-z), and H(2-z) is the Taylor series at 0 ===");
print("   [verified numerically in 03_values.py: H(z)+H(2-z)=2G to 40 digits]");
print("=== H is regular at z=0: H = sum a_n z^{n+1} in Q[[z]] ===");
print("=== conclusion: the singular set of H on P^1 is {1, infinity} ===");
print();
print("=== 2-adic radius of H in the two coordinates ===");
print("   z-coordinate:    v_2(a_n) = -1 exactly  =>  slope 0, R_2 = 1");
{
my(bad = 0);
for(n = 0, 300, if(valuation(2^(n+1)*amom(n), 2) != n, bad++));
print("   zeta-coordinate (z = 2 zeta^2): coefficient of zeta^{2n+2} is 2^{n+1} a_n,");
print("      v_2 = n exactly?  exceptions for n<=300: ", bad, "   => slope n/(2n+2) -> 1/2, R_2 = sqrt(2)");
}
quit;
