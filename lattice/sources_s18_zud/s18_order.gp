read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/s18_setup.gp");
B=vector(30); B[1]=0; B[2]=1;
A=vector(30); A[1]=1; A[2]=6;
for(n=1,28, B[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*B[n+1] - 12*n*(16*n^2-1)*B[n])/(n+1)^3);
for(n=1,28, A[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*A[n+1] - 12*n*(16*n^2-1)*A[n])/(n+1)^3);
print("A: ",vector(10,i,A[i]));
print("B: ",vector(10,i,B[i]));
print("F - sum A_n t^n = ", (F18 - sum(n=0,20,A[n+1]*T18^n))+O(x^21));
for(k=1,3, print("k=",k,": [t^n] F*theta^-k Phi = ", tcoeffs(F18*invth(PHI18,k),9)));
