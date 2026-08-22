\\ ================================================================
\\ The level-12 zeta(7) PARENT Apery system, built by the ZETA5_TWO_ROW recipe.
\\   source  Phi12 = C_{7,2,3}E_8, c=(1,-572,11583,-36608,46332,-20736) on d|12,
\\           L(Phi12,7) = (209/1728) zeta(7)
\\   host    x = eta(4t)^2 eta(12t)^2/(eta(t)^2 eta(3t)^2),  x|W_12 = 1/(16x)
\\           Fricke-invariant w = x/(1+16x^2); normalised t = w/(1+2w) = x/(16x^2+2x+1)
\\   rows    F = Phi12/(Dt) (weight 6), A(t) = F, B(t) = F * D^{-7}Phi12
\\ ================================================================
default(parisize,8000000000); default(realprecision,80);
N = 260;
et(k) = eta(q^k + O(q^(N+2)));
d12=[1,2,3,4,6,12]; c12=[1,-572,11583,-36608,46332,-20736];
{g = vector(N, m, sum(i=1,6, if(m%d12[i]==0, c12[i]*sigma(m\d12[i],7), 0)));}
Phi = sum(m=1,N, g[m]*q^m) + O(q^(N+1));
Psi = sum(m=1,N, (g[m]/m^7)*q^m) + O(q^(N+1));
xi  = (209/1728)*zeta(7);
x = q*et(4)^2*et(12)^2/(et(1)^2*et(3)^2);
t = x/(16*x^2+2*x+1);
print("x = ", x + O(q^7));
print("t = ", t + O(q^9));
Dt = q*deriv(t,q);
F  = Phi/Dt;
Qt = serreverse(t);
A  = subst(F,q,Qt);
B  = A*subst(Psi,q,Qt);
M  = N-10;
An = vector(M+1,i,polcoeff(A,i-1));
Bn = vector(M+1,i,polcoeff(B,i-1));
print("A_0..A_8 = ", vector(9,i,An[i]));
print("B_0..B_5 = ", vector(6,i,Bn[i]));
print("A_n integral to n=",M,"? ", vecmax(vector(M+1,i,denominator(An[i])))==1);
print("log|A_n|/n at n=",M," = ", log(abs(1.0*An[M+1]))/M, "   log(7+4sqrt3) = ", log(7+4*sqrt(3.)));
print("B_n/A_n - xi = ", 1.0*Bn[M+1]/An[M+1]-xi);
print("|B_n/A_n - xi|^(1/n) = ", abs(1.0*Bn[M+1]/An[M+1]-xi)^(1.0/M), "   (7-4sqrt3 = ",7-4*sqrt(3.),")");
print("lindep(B_M/A_M, zeta(7)) = ", lindep([1.0*Bn[M+1]/An[M+1], zeta(7)]));
\\ linear form
{Ln = vector(M+1,i, 1.0*Bn[i]-xi*An[i]);}
print("log|B_n - xi A_n|/n at n=",M," = ", log(abs(Ln[M+1]))/M);
\\ denominators
{my(dn=vector(M+1)); dn[1]=1; for(n=1,M, dn[n+1]=lcm(dn[n],n));
 for(k=4,8, my(ok=1,bad=0); for(n=1,M, if(denominator(dn[n+1]^k*Bn[n+1])!=1, ok=0;bad=n;break));
   print("  d_n^",k," B_n integral? ", ok, if(ok,"",concat("  first failure n=",bad))));}
print("v_ell(B_ell), ell=11..41: ", vector(9,i,my(l=[11,13,17,19,23,29,31,37,41][i]); valuation(Bn[l+1],l)));
\q
