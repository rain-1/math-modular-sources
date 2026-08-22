\\ Weight-two spaces at level 5 and just above: where a second class could live.
print("dim M_2(Gamma_0(5)) = ", mfdim([5,2,1],4), "   dim S_2 = ", mfdim([5,2,1],1));
print("M_2(Gamma_1(5)) split by character  [order, char, dim M_2, dim S_2]: ", mfdim([5,2,0],4));
print("  => trivial char: dim 1 (E_{2,5});  quadratic char chi_5: dim 2, all Eisenstein.");
print("dim M_2(Gamma_0(5),chi_5) = ", mfdim([5,2,[znstar(5,1),znconreychar(znstar(5,1),4)]],4));
print("");
print("higher levels on the same t-line (trivial character), dim M_2 / dim S_2:");
for(i=1,5, print("  N = ", 5*i, ":  ", mfdim([5*i,2,1],4), " / ", mfdim([5*i,2,1],1)));
print("");
print("Gamma_0(5)+5 has ONE cusp, so every weight-0 modular unit with cuspidal divisor");
print("has divisor 0 and is constant: the family F -> F*U is empty on this host.");
print("On Gamma_0(5) (two cusps) the unit is u=(eta_5/eta_1)^6, but a weight-2 form has");
print("divisor of degree mu/6 = 1, so ord_oo(E_{2,5}*u^j) = j forces j = 0.");
