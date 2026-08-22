default(parisizemax, 12000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/eta2000.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/eta500.gp");
print("v_5(eta) = ", valuation(eta2000, 5));
print("eta agrees between the N=500 and N=2000 runs to 5-adic precision 5^",
      valuation(eta2000 - eta500, 5));
e = eta2000;
print("eta = ", concat(Vecrev(apply(z->Str(z), digits(e % 5^25, 5)))), "  (base 5, most significant first, mod 5^25)");
print("eta 5-adic digits c_0,c_1,...,c_24 (LSB first): ", Vecrev(digits(e % 5^25, 5)));
print("eta mod 5^25 (decimal) = ", e % 5^25);
print("2*eta mod 5^25 (decimal) = ", (2*e) % 5^25);
print("eta as a 5-adic number: ", e + O(5^25));
print("2*eta as a 5-adic number (= zeta_5(3)?): ", 2*e + O(5^25));
quit;
