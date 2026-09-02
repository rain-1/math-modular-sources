\\ 01_setup.gp -- what PARI can do at level 28, weight 5/2 and 29/2; cusps; Delta(4tau).
default(parisize, 800000000);
print("=== mfinit([28,5/2]) ===");
M = mfinit([28,5/2]);
print("dim M_{5/2}(28) = ", mfdim(M));
S = mfinit([28,5/2],1);
print("dim S_{5/2}(28) = ", mfdim(S));
K = mfkohnenbasis(S);
print("kohnen basis matrix size = ", matsize(K));
B = mfbasis(S);
print("num cusp basis = ", #B);
f = mflinear(S, K[,1]);
print("kohnen form coeffs 0..40: ", mfcoefs(f,40));
print("=== cusps of 28 ===");
print(mfcusps(28));
print("=== dim M_{29/2}(28) ===");
M2 = mfinit([28,29/2]);
print("dim = ", mfdim(M2));
quit;
