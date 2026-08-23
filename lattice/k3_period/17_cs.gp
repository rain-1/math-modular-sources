default(parisizemax, 8000000000);
default(realprecision, 160);
mf32 = mfinit([32,3,-8],0); gg = mfeigenbasis(mf32)[1];
LL = lfunmf(mf32, gg);
L1 = lfun(LL,1); L2 = lfun(LL,2);
print("L(g,1) = ", L1);
print("L(g,2) = ", L2);
print("L(g,2)/L(g,1) = ", L2/L1);
print("L(g,2)/(Pi*L(g,1)) = ", L2/(Pi*L1));
/* Chowla-Selberg period for discriminant -8 (class number 1, w=2) */
om = (gamma(1/8)*gamma(3/8)/(gamma(5/8)*gamma(7/8)))^(1/2);
print();
print("CS quantity om = (G(1/8)G(3/8)/(G(5/8)G(7/8)))^(1/2) = ", om);
/* Deligne: for weight 3 CM by Q(sqrt-2), L(g,2) ~ Pi^2 * (period)^2 ... probe */
{
foreach([1,2,3,4], e,
  print("L(g,2)/(Pi^", e, ") = ", L2/Pi^e));
}
print();
/* candidate basis lindep for L(g,2) */
bas = [L2, Pi^2*om^2, Pi*om^2, om^4, Pi^2, Pi^3, sqrt(2)*Pi^2, Pi^2/sqrt(2)];
print("lindep(L2, Pi^2 om^2, Pi om^2, om^4, ...) = ", lindep(bas));
print("lindep([L2, Pi^2*om^2]) = ", lindep([L2, Pi^2*om^2]));
print("lindep([L2, Pi*om^2])  = ", lindep([L2, Pi*om^2]));
print("lindep([L2, om^4])     = ", lindep([L2, om^4]));
print("L2/(Pi^2*om^2) = ", L2/(Pi^2*om^2));
print("L2/(Pi*om^2) = ", L2/(Pi*om^2));
print("L2/om^4 = ", L2/om^4);
print();
print("L(g,1) vs Pi*om^2 : ", L1/(Pi*om^2), "  lindep ", lindep([L1, Pi*om^2]));
print("L(g,1)/L(g,2) * Pi = ", Pi*L1/L2);
quit;
