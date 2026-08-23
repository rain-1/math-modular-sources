/* (b) the inhomogeneous ODE for Phi and the collapsed H-equation.
   Names avoided: Phi -> phser, H -> hser, S -> uni. */
default(parisizemax, 4*10^9);
NN = 300;
amom(n) = sum(k = 0, n, binomial(n,k)/(2*k+1)) / (2^(n+1)*(n+1));
av = vector(NN+1); for(n = 0, NN, av[n+1] = amom(n));
aa(n) = if(n < 0, 0, av[n+1]);

print("=== (b1) theta-operator form: [(t+1)(2t+1) - z(t+1)(3t+2) + z^2(t+1)^2] Phi = 1/2 ===");
{
my(bad = 0);
/* coefficient of z^m in the LHS */
for(m = 0, NN-2,
  my(c = (m+1)*(2*m+1)*aa(m) - if(m>=1, m*(3*m-1)*aa(m-1), 0) + if(m>=2, (m-1)^2*aa(m-2), 0));
  my(want = if(m == 0, 1/2, 0));
  if(c != want, bad++; if(bad<4, print("FAIL m=",m," got ",c," want ",want)))
);
print("theta-form failures m<=", NN-2, ": ", bad);
}

print("=== (b2) collapsed form z(z-2)H'' + (z-1)H' = 1/(2(z-1)), H = z*Phi ===");
print("    coefficientwise this is  m^2 a_{m-1} - (m+1)(2m+1) a_m = -1/2  for all m>=0");
{
my(bad = 0);
for(m = 0, NN,
  my(c = m^2*aa(m-1) - (m+1)*(2*m+1)*aa(m));
  if(c != -1/2, bad++; if(bad<4, print("FAIL m=",m," got ",c)))
);
print("collapsed-form failures m<=", NN, ": ", bad);
}

print("=== (b3) the two-term inhomogeneous recurrence IMPLIES the note's three-term one ===");
print("    (n+3)(2n+5)a_{n+2} - (n+2)(3n+5)a_{n+1} + (n+1)^2 a_n");
print("    = [ (n+3)(2n+5)a_{n+2} - (n+2)^2 a_{n+1} ] - (n+2)[ (n+2)... ]  -- checked numerically in 01");

print("=== (b4) homogeneous pure part: 1 and uni(z) = 2*asin(sqrt(z/2)) ===");
print("    substitute z = 2*w^2 so that uni = 2*asin(w) is a power series in w;");
print("    then d/dz = (1/(4w)) d/dw.");
{
my(w = 'w + O('w^30));
my(uni = 2*asin(w));
my(d1 = deriv(uni,'w)/(4*'w));
my(d2 = deriv(d1,'w)/(4*'w));
my(zz = 2*'w^2);
print("    residual of z(z-2)Y''+(z-1)Y' on Y = uni: ", zz*(zz-2)*d2 + (zz-1)*d1);
print("    (Y = 1 is trivially a solution.)");
}
quit;
