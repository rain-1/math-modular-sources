default(parisize,"4G");
read("lib.gp");
NT=30;
M = mirror((n+1)^3, (2*n+1)*(17*n^2+17*n+5), -1*n^3, NT);
print("gamma (Apery 17,5,1):");
print("  t(q) = ", M[1]+O('t^9));
print("  F(q) = ", M[2]+O('t^9));
print("  y0(t)= ", M[3]+O('t^7));
E(k,M)=eta(q^k+O(q^(M+2)))
print("  eta check t: ", (q*prod(j=1,NT+2,((1-q^j+O(q^(NT+2)))*(1-q^(6*j)+O(q^(NT+2))))/((1-q^(2*j)+O(q^(NT+2)))*(1-q^(3*j)+O(q^(NT+2)))))^12) + O(q^9));
quit;
