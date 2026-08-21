default(parisizemax,2000000000);
MM=260;
default(seriesprecision,MM+4);
et(k)=eta(q^k+O(q^(MM+4)));
F=et(2)^7*et(3)^7/(et(1)^5*et(6)^5); t=q*(et(1)*et(6)/(et(2)*et(3)))^12;
f=sqrt(F); Psi=f^3*t/4;
\\ check Psi = (1/4) (eta1 eta6)^{9/2} (eta2 eta3)^{-3/2}
Psi2=(q*et(1)*et(6))^(9/2)/(et(2)*et(3))^(3/2)/4;  \\ q-power: (9/2*(1+6) - 3/2*(2+3))/24 = (31.5-7.5)/24 = 1 -> q^1 factor
print("Psi identity: ",Psi-Psi2==0);
print("4^m psi(m) integral? ",vector(MM-2,m,denominator(4^m*polcoeff(Psi,m)))==vector(MM-2,m,1));
Th=sum(m=1,MM-2,polcoeff(Psi,m)/m^2*q^m);
comp=f*Th; tt=serreverse(t); comp_t=subst(comp,q,tt);
b_from_identity=vector(MM-4,n,4^(n-1)*polcoeff(comp_t,n-1));
N=MM-4; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=10;B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=((136*n^2+68*n+10)*A[n+1]-4*(2*n-1)^2*A[n])/(n+1)^2;B[n+2]=((136*n^2+68*n+10)*B[n+1]-4*(2*n-1)^2*B[n])/(n+1)^2);
print("b_n from identity == recurrence companion for n<",N,": ",vector(N,i,b_from_identity[i])==vector(N,i,B[i]));
\\ a_n from f as well
ft=subst(f,q,tt); print("a_n check: ",vector(N,i,4^(i-1)*polcoeff(ft,i-1))==vector(N,i,A[i]));
\q
