/* 03_ode.gp -- L = x(1-11x-x^2) y'' + (1-22x-3x^2) y' - (3+x) y.
   [x^(n-1)] L(y) = n^2 y_n - (11n^2-11n+3) y_{n-1} - (n-1)^2 y_{n-2}.
   Compute the RHS series for A, BD, B3, B4 exactly and identify it.        */
default(parisizemax, 20000000000);
read("build.txt");
Lc(bb, n) = my(b0=bb[n+1], b1=if(n>=1, bb[n], 0), b2=if(n>=2, bb[n-1], 0)); n^2*b0 - (11*n^2-11*n+3)*b1 - (n-1)^2*b2;
show(bb, nm) = print(nm, ":  [x^0..x^11] L = ", vector(12, i, Lc(bb, i)));
show(A, "L A ");
show(BD, "L BD");
show(B3, "L B3");
show(B4, "L B4");
/* check homogeneity of A to high n */
mx = 0; for(n=1, NA, if(Lc(A,n)!=0, mx=1)); print("L A identically 0 up to n=", NA, "? ", mx==0);
/* series of RHS for each */
rD = vector(NA, n, Lc(BD, n));
r3 = vector(NA, n, Lc(B3, n));
r4 = vector(NA, n, Lc(B4, n));
print("");
print("rD = ", vector(14,i,rD[i]));
print("r3 = ", vector(14,i,r3[i]));
print("r4 = ", vector(14,i,r4[i]));
/* denominators of 1/(1-11x-x^2)-type expansions: test rD = alpha + beta*x/(1-11x-x^2) etc.
   Fit RHS = (p0+p1 x+p2 x^2 + ...)/(1-11x-x^2) with small numerator.        */
fit(r, nm, deg) = my(S = sum(i=1,NA, r[i]*x^(i-1)) + O(x^NA)); my(N = S*(1-11*x-x^2)); print(nm, " * (1-11x-x^2) = ", vector(deg+4, i, polcoeff(N,i-1)));
fit(rD, "rD", 6);
fit(r3, "r3", 6);
fit(r4, "r4", 6);
quit;
