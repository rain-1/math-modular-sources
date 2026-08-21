PS=100; read("cover.gp");
chi(n)=kronecker(-3,n);
S = sum(m=1,PREC-1, sumdiv(m,d, chi(m/d)*d^2)*q^m) + O(q^PREC);
Psi0 = S - 4*Vd(S,2);
ThPsi = Dqinv(Psi0,2);
ThC = Dqinv(PhiC,2);
rev = serreverse(tC);              \\ q as series in t
Aser = subst(FC,q,rev);            \\ = sum a^C_n t^n
Bser = subst(FC*ThC,q,rev);        \\ = sum b^C_n t^n
Pser = subst(FC*ThPsi,q,rev);      \\ companion of Psi0 on C's curve
print("a_n check (should be 1,3,9,...): ",vector(6,i,polcoeff(Aser,i-1)));
print("b_n check (0,1,...): ",vector(6,i,polcoeff(Bser,i-1)));
print("psi_n : ",vector(6,i,polcoeff(Pser,i-1)));
print("n  v3(psi_n)  v3(a_n)  v3(psi_n/a_n)");
for(k=1,9, my(n=10*k); print(n,"  ",valuation(polcoeff(Pser,n),3),"  ",valuation(polcoeff(Aser,n),3),"  ",valuation(polcoeff(Pser,n)/polcoeff(Aser,n),3)));
