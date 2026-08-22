default(parisizemax, 4000000000);
default(realprecision, 60);
NK = 150;
LCT = vector(NK+2); LCT[1]=1; for(n=1,NK+1, LCT[n+1]=lcm(LCT[n],n));
badat(b,k,N) = my(bad=-1); for(n=1,N, if(denominator(LCT[n+1]^k*b[n+1])!=1, bad=n; break)); bad
sharpk(b,N) = my(kk=-1); for(k=0,8, if(badat(b,k,N)<0, kk=k; break)); kk
P0(n)=(n+2)^3
P1(n)=22*n^3+165*n^2+419*n+360
P2(n)=119*n^3+1071*n^2+3222*n+3240
P3(n)=-(22*n^3+231*n^2+815*n+966)
P4(n)=(n+4)^3
R1(n)=2*(2*n-1)^2*(11*n^2-11*n+5)
R2(n)=500*(2*n-3)*(n-1)^2*(2*n-1)
\\ =================== 12. Sym^2 of Zagier D ==============================
N=260;
a = vector(N+1); a[1]=1; a[2]=3; for(n=1,N-1, a[n+2]=((11*n^2+11*n+3)*a[n+1]+n^2*a[n])/(n+1)^2);
AA = vector(N+1); for(n=0,N, AA[n+1]=sum(i=0,n, a[i+1]*a[n-i+1]));
print("=== 12. Sym^2 of Zagier D");
print("  Zagier D a_n : ", vector(8,i,a[i]));
print("  Sym^2  A_n   : ", vector(8,i,AA[i]));
print("  A_n integral for n<=",N," : ", if(sum(n=0,N,denominator(AA[n+1])!=1)==0,"yes","NO"));
chk = sum(t=0,N-5, abs(P0(t)*AA[t+1]+P1(t)*AA[t+2]+P2(t)*AA[t+3]+P3(t)*AA[t+4]+P4(t)*AA[t+5]));
print("  CDT_FINDER 5-term recurrence residual sum over 0<=n<=",N-5," : ", chk);
cp12 = x^4-22*x^3+119*x^2+22*x+1;
print("  charpoly = ", cp12, "   factor over Q : ", factor(cp12));
print("  roots : ", polroots(cp12));
phi=(1+sqrt(5))/2; print("  phi^5 = ", phi^5, "  -phi^-5 = ", -phi^-5);
BB = vector(N+1); BB[1]=0; BB[2]=1; BB[3]=0; BB[4]=0;
for(t=0,N-5, BB[t+5] = -(P0(t)*BB[t+1]+P1(t)*BB[t+2]+P2(t)*BB[t+3]+P3(t)*BB[t+4])/P4(t));
k12 = sharpk(BB,NK);
print("  companion (b0,b1,b2,b3)=(0,1,0,0): sharp k (n<=",NK,") = ", k12, "   first fail at k-1: n=", badat(BB,k12-1,NK));
print("  B_n/A_n at n=200 : ", BB[201]*1.0/AA[201]);
l2 = 0.0901699437494742410229341718281905886015458990288143106772431;
print("  second-largest DISTINCT |root| = ", l2, "  rho = ", 1/l2, "  margin = ", -log(l2)-k12);
print("  with-multiplicity reading lambda_2 = phi^5 : margin = ", -log(phi^5)-k12);
print("");
\\ =================== 13. AESZ 184 =======================================
print("=== 13. AESZ 184");
N2=200;
u = vector(N2+1); u[1]=1; u[2]=10; for(n=2,N2, u[n+1]=(R1(n)*u[n]-R2(n)*u[n-1])/n^4);
print("  u_0..u_6 : ", vector(7,i,u[i]));
print("  u_n integral for n<=",N2," : ", if(sum(n=0,N2,denominator(u[n+1])!=1)==0,"yes","NO"));
b = vector(N2+1); b[1]=0; b[2]=1; for(n=2,N2, b[n+1]=(R1(n)*b[n]-R2(n)*b[n-1])/n^4);
k13 = sharpk(b,NK);
print("  companion b0=0,b1=1: sharp k (n<=",NK,") = ", k13, "  first fail at k-1: n=", badat(b,k13-1,NK));
cp13 = x^2-88*x+2000;
print("  charpoly = ", cp13, "  factor: ", factor(cp13), "  disc = ", poldisc(cp13));
print("  roots : ", polroots(cp13), "   |root| = sqrt(2000) = ", sqrt(2000));
print("  margin = -log|l2| - k = ", -log(sqrt(2000))-k13, "   (complex conjugate pair: no fold)");
print("  b_n/u_n at n=200 : ", b[201]*1.0/u[201]);
