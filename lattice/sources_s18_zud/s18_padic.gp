read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
PR = 220; N = 340;
A=vector(N+2); B=vector(N+2);
A[1]=1; A[2]=6; B[1]=0; B[2]=1;
{ for(n=1,N,
    A[n+2]=(2*(2*n+1)*(7*n^2+7*n+3)*A[n+1]-12*n*(16*n^2-1)*A[n])/(n+1)^3;
    B[n+2]=(2*(2*n+1)*(7*n^2+7*n+3)*B[n+1]-12*n*(16*n^2-1)*B[n])/(n+1)^3); }
xi = B[N+1]/A[N+1];  xim = B[N]/A[N];
print("N=",N-1,"  Cauchy v_3(xi_N-xi_{N-1}) = ", valuation(xi-xim,3));
print("v_3(A_N)=",valuation(A[N+1],3),"  v_3(B_N)=",valuation(B[N+1],3));
z2 = Lp(3,triv,2,PR);  z3 = Lp(3,triv,3,PR);
print("v_3(xi - zeta_3(2)/2) = ", valuation(xi - z2/2, 3));
print("v_3(xi - zeta_3(3)/2) = ", valuation(xi - z3/2, 3));
print("v_3(xi + zeta_3(2)/2) = ", valuation(xi + z2/2, 3));
quit;
