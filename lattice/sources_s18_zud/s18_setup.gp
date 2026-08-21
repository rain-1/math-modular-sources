/* Cooper s18: t,F,Phi = F*theta_q t.  Shared setup. */
NQ = 46;
ee(d) = eta(x^d + O(x^NQ));
E2s(d) = 1 - 24*sum(n=1,(NQ-1)\d, sigma(n)*x^(d*n)) + O(x^NQ);
U18 = ee(3)^4*ee(6)^4/(ee(1)^2*ee(2)^2*ee(9)^2*ee(18)^2);   /* = q*x_18 */
T18 = x*U18/(U18+3*x)^2;
F18 = (1/4)*(18*E2s(18) - 9*E2s(9) - 12*E2s(6) + 6*E2s(3) + 2*E2s(2) - E2s(1));
PHI18 = F18*x*deriv(T18);
M = NQ-3;
invth(f,k) = sum(m=1,M, polcoeff(f,m)/m^k*x^m) + O(x^(M+1));
tcoeffs(Y,nmax) = {my(co=vector(nmax+1),R=Y); for(n=0,nmax, co[n+1]=polcoeff(R,0); R=R-co[n+1]; if(n<nmax, R=R/T18)); co}
