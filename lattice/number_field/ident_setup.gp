default(realprecision, 120);
\\ Dirichlet L-values for Kronecker characters
LK(D, s) = { my(q=abs(D), t=0); for(j=1,q-1, t += kronecker(D,j)*zetahurwitz(s, j/q)); t/q^s; }
print("L(chi_-7,1) = ", LK(-7,1), "   pi/sqrt7 = ", Pi/sqrt(7));
print("L(chi_-7,2) = ", LK(-7,2));
print("L(chi_-3,1) = ", LK(-3,1), "   pi/(3sqrt3) = ", Pi/(3*sqrt(3)));
print("L(chi_-3,2) = ", LK(-3,2));
\\ newforms
mf7 = mfinit([7,3,-7],0); print("dim S_3^new(7,chi-7) = ", mfdim(mf7));
E7 = mfeigenbasis(mf7); print("#eigen = ", #E7);
f7 = E7[1]; print("f7 coeffs: ", mfcoefs(f7, 12));
L7 = lfunmf(mf7, f7);
print("L(f7,1) = ", lfun(L7,1));
print("L(f7,2) = ", lfun(L7,2));
mf9 = mfinit([9,3,-3],0); print("dim S_3^new(9,chi-3) = ", mfdim(mf9));
E9 = mfeigenbasis(mf9); print("#eigen = ", #E9);
f9 = E9[1]; print("f9 coeffs: ", mfcoefs(f9, 12));
L9 = lfunmf(mf9, f9);
print("L(f9,1) = ", lfun(L9,1));
print("L(f9,2) = ", lfun(L9,2));
