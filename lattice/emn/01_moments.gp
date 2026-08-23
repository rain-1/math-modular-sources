/* EMN verification: claim (a) -- moments a_n and their recurrence.
   a_n = \iint_Delta (x^2+y^2)^n dx dy, Delta = {x,y>=0, x+y<=1}.
   Note's closed form: a_n = 1/(2^{n+1}(n+1)) * sum_k binom(n,k)/(2k+1).
   Recurrence: (n+3)(2n+5) a_{n+2} - (n+2)(3n+5) a_{n+1} + (n+1)^2 a_n = 0.
   Avoid builtin names: use amom, hh, ss.  */
default(parisizemax, 4*10^9);

/* direct exact double integral of (x^2+y^2)^n over the simplex */
amom_direct(n) = {
  my(t = 0);
  /* (x^2+y^2)^n = sum_j binom(n,j) x^{2j} y^{2n-2j} */
  for(j = 0, n,
    /* \int_0^1 \int_0^{1-x} x^{2j} y^{2n-2j} dy dx
       = \int_0^1 x^{2j} (1-x)^{2n-2j+1}/(2n-2j+1) dx
       = B(2j+1, 2n-2j+2)/(2n-2j+1) */
    t += binomial(n,j) * (2*j)! * (2*n-2*j+1)! / (2*n+2)! / (2*n-2*j+1)
  );
  t;
};

amom_closed(n) = {
  my(s = sum(k = 0, n, binomial(n,k)/(2*k+1)));
  s / (2^(n+1)*(n+1));
};

print("=== (a1) closed form vs direct integral, n<=60 ===");
{
my(bad = 0);
for(n = 0, 60,
  if(amom_direct(n) != amom_closed(n), bad++; print("MISMATCH n=", n))
);
print("closed-form mismatches n<=60: ", bad);
}

print("=== (a2) recurrence check, n<=200 ===");
{
my(NN = 205, av = vector(NN+1), bad = 0, worst = -1);
for(n = 0, NN, av[n+1] = amom_closed(n));
for(n = 0, NN-2,
  my(r = (n+3)*(2*n+5)*av[n+3] - (n+2)*(3*n+5)*av[n+2] + (n+1)^2*av[n+1]);
  if(r != 0, bad++; if(worst<0, worst=n); print("REC FAIL n=", n, " residual=", r))
);
print("recurrence failures for 0<=n<=", NN-2, ": ", bad);
}

print("=== (a3) numeric a_n asymptotics ===");
{
for(n = 1, 6, print("a_", n, " = ", amom_closed(n)));
print("a_50 ~ ", 1.0*amom_closed(50));
print("a_100 ~ ", 1.0*amom_closed(100));
print("ratio a_101/a_100 = ", 1.0*amom_closed(101)/amom_closed(100));
print("ratio a_201/a_200 = ", 1.0*amom_closed(201)/amom_closed(200));
}
quit;
