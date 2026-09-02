\\ Fricke eigenvalues of W_N on S_4(Gamma_0(N)) for the Fricke-family host levels.
default(parisizemax, 8000000000);
NC = 40;
doN(N) = my(mf, d, B, ai, M, Kp, Km, g, v); mf = mfinit([N,4],1); d = mfdim(mf); print("### N=", N, "  dimS4 = ", d); if(d==0, return(0)); B = mfbasis(mf); ai = mfatkininit(mf, N); M = matrix(d,d); for(j=1,d, v = mftobasis(mf, mfatkin(ai, B[j])); for(t=1,d, M[t,j] = v[t])); print("  W_N matrix = ", M); print("  charpoly = ", charpoly(M)); Kp = matker(M-1); Km = matker(M+1); print("  dim(+1) = ", matsize(Kp)[2], "   dim(-1) = ", matsize(Km)[2]); for(j=1,matsize(Km)[2], g = mflinear(mf, Km[,j]); print("  MINUS[",j,"] = ", mfcoefs(g,NC))); for(j=1,matsize(Kp)[2], g = mflinear(mf, Kp[,j]); print("  PLUS [",j,"] = ", mfcoefs(g,NC))); 1;
lev = [5,6,7,8,9,10,12,18];
for(i=1,#lev, doN(lev[i]));
quit;
