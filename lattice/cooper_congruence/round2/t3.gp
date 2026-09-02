read("wt2.gp");
default(realprecision,50);
t = 0.11 + 0.29*I;
\\ gamma in Gamma_0(7):  [3,1;7,...]?  need ad-bc=1 with 7|c: [a,b;c,d]=[1,0;7,1]
g1 = [1,0;7,1];
tt = (g1[1,1]*t+g1[1,2])/(g1[2,1]*t+g1[2,2]);
A = fhat7(t); B = fhat7(tt);
print("f(t)  = ", A[1]);
print("f(gt)*(7t+1)^2 = ", B[1]*(7*t+1)^2, "   (should equal f(t): weight -2)");
print("fhat(t)  = ", A[2]);
print("fhat(gt) = ", B[2], "   (should equal fhat(t): weight 0 invariant)");
g2 = [8,1;7,1];
tt2 = (8*t+1)/(7*t+1);
C = fhat7(tt2);
print("fhat(g2 t) = ", C[2]);
quit;
