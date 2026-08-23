/* Task 2(3): the 3-adic pair.
   EMN 3-adic class  L_3(2,chi_12)  vs  the conductor-6 decayer's 3-adic class zeta_3(2). */
default(parisizemax, 8*10^9);
PREC = 80;
chr(a) = my(r = a % 12); if(r==1 || r==11, 1, if(r==5 || r==7, -1, 0));
tri(a) = if(a % 3 == 0, 0, 1);            /* trivial character */
teich(a) = if(a % 3 == 1, 1, if(a % 3 == 2, -1, 0));

lkl(sv, FF, JMAX, which) = my(tot = 0); for(a = 1, FF, my(c = if(which==12, chr(a), tri(a))); if(a % 3 != 0 && c != 0, my(br = a/teich(a)); my(inner = sum(j = 0, JMAX, binomial(1-sv,j)*(FF/a)^j*bernfrac(j))); tot += c*(br + O(3^PREC))^(1-sv)*inner)); tot/(FF*(sv-1));

L3chi12 = lkl(2, 12, 130, 12);
z3 = lkl(2, 3, 130, 1);                   /* zeta_3(2) = L_3(2, trivial) */
print("=== the two 3-adic classes ===");
print("  L_3(2,chi_12) = ", lift(L3chi12 + O(3^18)), " + O(3^18)");
print("  zeta_3(2)     = ", lift(z3 + O(3^18)), " + O(3^18)");
print();
print("  sanity: zeta_3(-1) should be -(1-3^0)B_2/2 = 0  (trivial character, Euler factor kills it)");
print("    zeta_3(0)  = ", lift(lkl(0,3,130,1) + O(3^12)), "   (predicted -(1-1/3)B_1... )");
print("    zeta_3(-1) = ", lift(lkl(-1,3,130,1) + O(3^12)), "   (predicted -(1-3)B_2/2 = 1/6 = ", lift(1/6 + O(3^12)), ")");
print("    zeta_3(-3) = ", lift(lkl(-3,3,130,1) + O(3^12)), "   (predicted -(1-3^3)B_4/4 = ", lift(-(1-27)*(-1/30)/4 + O(3^12)), ")");
print();
print("=== are they 3-adically Q-proportional? ===");
{
my(hit = 0);
for(p = -40, 40, for(q = 1, 40, if(gcd(abs(p),q)==1 && p != 0,
   my(d = L3chi12 - (p/q)*z3); if(valuation(d,3) > 40, print("  L_3(2,chi_12) = (",p,"/",q,") zeta_3(2)  to v_3=", valuation(d,3)); hit++))));
if(!hit, print("  no rational p/q with |p|,q <= 40 gives v_3 > 40  =>  NOT proportional"));
}
print("  lindep on [L_3(2,chi_12), zeta_3(2), 1] to 3-adic precision 78:");
print("   ", lindep([L3chi12 + O(3^78), z3 + O(3^78), 1 + O(3^78)]));
print("  (a huge-height relation means: 3-adically independent at this precision)");
print();
print("=== the EMN 3-adic row  s_N = sum_{n<=N} a_n (3/2)^n ===");
amom(n) = sum(k = 0, n, binomial(n,k)/(2*k+1))/(2^(n+1)*(n+1));
{
my(s = 0, sv = vector(121));
for(n = 0, 120, s += amom(n)*(3/2)^n; sv[n+1] = s);
print("  N   v_3(s_N - L_3)   log|num s_N|/N   log den(s_N)/N   v_2(den s_N)/N");
foreach([20,40,60,80,100,120], nn,
  my(x = sv[nn+1], d = denominator(x));
  print("  ", nn, "   ", valuation(x - L3chi12, 3), "        ",
        log(abs(numerator(x)*1.0))/nn, "   ", log(1.0*d)/nn, "   ", valuation(d,2)/nn));
}
print();
print("=== the EMN archimedean row  t_N = sum_{n<=N} a_n 2^{-n}  ->  (2/3) G ===");
{
my(s = 0);
my(gg = 0.9159655941772190150546035149323841107741);
for(n = 0, 200, s += amom(n)/2^n; if(n==20 || n==50 || n==100 || n==150 || n==200,
  my(d = denominator(s));
  print("  N=", n, "  -log|2G/3 - t_N|/N = ", -log(abs(2*gg/3 - s*1.0))/n,
        "   log den(t_N)/N = ", log(1.0*d)/n, "   v_2(den)/N = ", 1.0*valuation(d,2)/n)));
}
quit;
