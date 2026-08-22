default(parisizemax, 4000000000);
default(realprecision, 1200);
\\ --- AESZ 207 : two-term linear form alpha*A_n - B_n ------------------
Q0(t)=t^4
Q1(t)=-2^4*(1072*t^4 - 17824*t^3 - 10888*t^2 - 1976*t - 145)
Q2(t)=-2^17*(51088*t^4 + 116368*t^3 - 45264*t^2 - 14228*t - 1397)
Q3(t)=2^28*13*(73104*t^4 + 1536*t^3 - 488*t^2 + 384*t + 97)
Q4(t)=-2^44*13^2*(2*t+1)^4
N=460;
A = vector(N+1); A[1]=1.0;
for(n=1,N, A[n+1] = -( if(n>=1,Q1(n-1)*A[n],0) + if(n>=2,Q2(n-2)*A[n-1],0) + if(n>=3,Q3(n-3)*A[n-2],0) + if(n>=4,Q4(n-4)*A[n-3],0) )/Q0(n));
B = vector(N+1); B[1]=0.0; B[2]=1.0; B[3]=0.0; B[4]=0.0;
for(n=4,N, B[n+1] = -( Q1(n-1)*B[n] + Q2(n-2)*B[n-1] + Q3(n-3)*B[n-2] + Q4(n-4)*B[n-3] )/Q0(n));
al = B[N+1]/A[N+1];
r = vector(N+1, i, al*A[i]-B[i]);
default(realprecision,25);
print("AESZ 207: alpha = lim B_n/A_n = ", al);
print("  two-term form r_n = alpha A_n - B_n, log|r_n/r_{n-1}| at n=100,200,300,380 :");
print("   ", log(abs(r[101]/r[100])), " ", log(abs(r[201]/r[200])), " ", log(abs(r[301]/r[300])), " ", log(abs(r[381]/r[380])));
print("  log 53248 = ", log(53248.0), "   log 187.389 = ", log(187.3892067201468), "   log 89531.389 = ", log(89531.3892067201468));
default(realprecision,1200);
print("");
\\ --- Sym^2 Zagier D : two-term linear form -----------------------------
P0(n)=(n+2)^3
P1(n)=22*n^3+165*n^2+419*n+360
P2(n)=119*n^3+1071*n^2+3222*n+3240
P3(n)=-(22*n^3+231*n^2+815*n+966)
P4(n)=(n+4)^3
M=460;
a = vector(M+1); a[1]=1.0; a[2]=3.0; for(n=1,M-1, a[n+2]=((11*n^2+11*n+3)*a[n+1]+n^2*a[n])/(n+1)^2);
AA = vector(M+1); for(n=0,M, AA[n+1]=sum(i=0,n, a[i+1]*a[n-i+1]));
BB = vector(M+1); BB[1]=0.0; BB[2]=1.0; BB[3]=0.0; BB[4]=0.0;
for(t=0,M-5, BB[t+5] = -(P0(t)*BB[t+1]+P1(t)*BB[t+2]+P2(t)*BB[t+3]+P3(t)*BB[t+4])/P4(t));
be = BB[M+1]/AA[M+1];
rr = vector(M+1, i, be*AA[i]-BB[i]);
default(realprecision,25);
print("Sym^2 Zagier D: beta = lim B_n/A_n = ", be);
print("  two-term form log|r_n/r_{n-1}| at n=100,200,300,380 :");
print("   ", log(abs(rr[101]/rr[100])), " ", log(abs(rr[201]/rr[200])), " ", log(abs(rr[301]/rr[300])), " ", log(abs(rr[381]/rr[380])));
print("  log(phi^5) = ", log(11.09016994374947424102293), "   log(phi^-5)= ", -log(11.09016994374947424102293));
print("");
\\ --- weight-1 Zagier D two-term form (control) --------------------------
b1 = vector(M+1); b1[1]=0.0; b1[2]=1.0; for(n=1,M-1, b1[n+2]=((11*n^2+11*n+3)*b1[n+1]+n^2*b1[n])/(n+1)^2);
g = b1[M+1]/a[M+1];
r1 = vector(M+1, i, g*a[i]-b1[i]);
print("Zagier D (weight 1) control: gamma = lim b_n/a_n = ", g, "   zeta(2)/5 = ", zeta(2)/5);
print("  log|r_n/r_{n-1}| at n=100,200,300,380 : ", log(abs(r1[101]/r1[100])), " ", log(abs(r1[201]/r1[200])), " ", log(abs(r1[301]/r1[300])), " ", log(abs(r1[381]/r1[380])));
