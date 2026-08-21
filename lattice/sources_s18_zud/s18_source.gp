\\ Cooper s18: modular parametrisation, source Phi = F * theta_q t, Eisenstein decomposition
default(parisizemax,4000000000);
N = 60;  \\ q-precision
ee(d) = eta(x^d + O(x^N));           \\ prod (1-q^{dn})
E2s(d) = { my(v); v = 1 + O(x^N); v = 1 - 24*sum(n=1,N-1, sigma(n)*x^(d*n)); v; }  \\ E_2(d tau)

\\ x_18 = q^{-1} * e3^4 e6^4 / (e1^2 e2^2 e9^2 e18^2)
X = ee(3)^4*ee(6)^4/(ee(1)^2*ee(2)^2*ee(9)^2*ee(18)^2);
\\ work with u = q*x_18  (a power series with constant term 1)
U = X;    \\ this is q*x_18
\\ t = x/(x+3)^2 = (U/q)/((U+3q)/q)^2 = q*U/(U+3q)^2
T = x*U/(U+3*x)^2;
print("t_18 = ", T + O(x^12));

F = (1/4)*(18*E2s(18) - 9*E2s(9) - 12*E2s(6) + 6*E2s(3) + 2*E2s(2) - E2s(1));
print("F_18 = ", F + O(x^12));

\\ check F = sum s18(n) t^n
s18 = vector(20); s18[1]=1;  \\ s18[k] = A_{k-1}
A=vector(21); A[1]=1; A[2]=6;
for(n=1,19, A[n+2] = (2*(2*n+1)*(7*n^2+7*n+3)*A[n+1] - 12*n*(16*n^2-1)*A[n])/(n+1)^3);
print("A_n = ", A);
chk = sum(n=0,15, A[n+1]*T^n) - F;
print("F - sum A_n t^n = ", chk + O(x^16));

\\ Wronskian check: theta_q log t = F * sqrt(1-28t+192t^2)
Dq(f) = x*deriv(f);
K = Dq(T)/T;
print("K - F*sqrt(P) = ", (K - F*sqrt(1-28*T+192*T^2)) + O(x^16));

\\ The source
PHI = F*Dq(T);
print("Phi = ", PHI + O(x^30));
print("Phi coeffs: ", Vec(PHI + O(x^40)));
