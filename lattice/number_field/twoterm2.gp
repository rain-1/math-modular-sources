default(parisizemax, 4000000000);
default(realprecision, 2600);
P0(n)=(n+2)^3
P1(n)=22*n^3+165*n^2+419*n+360
P2(n)=119*n^3+1071*n^2+3222*n+3240
P3(n)=-(22*n^3+231*n^2+815*n+966)
P4(n)=(n+4)^3
M=420;
a = vector(M+1); a[1]=1.0; a[2]=3.0; for(n=1,M-1, a[n+2]=((11*n^2+11*n+3)*a[n+1]+n^2*a[n])/(n+1)^2);
AA = vector(M+1); for(n=0,M, AA[n+1]=sum(i=0,n, a[i+1]*a[n-i+1]));
BB = vector(M+1); BB[1]=0.0; BB[2]=1.0; BB[3]=0.0; BB[4]=0.0;
for(t=0,M-4, BB[t+5] = -(P0(t)*BB[t+1]+P1(t)*BB[t+2]+P2(t)*BB[t+3]+P3(t)*BB[t+4])/P4(t));
be = BB[M+1]/AA[M+1];
rr = vector(M+1, i, be*AA[i]-BB[i]);
default(realprecision,25);
print("Sym^2 Zagier D: beta = lim B_n/A_n = ", be);
print("  two-term form log|r_n/r_{n-1}| at n=100,200,300,400 : ", log(abs(rr[101]/rr[100])), " ", log(abs(rr[201]/rr[200])), " ", log(abs(rr[301]/rr[300])), " ", log(abs(rr[401]/rr[400])));
print("  log(phi^5) = ", log(11.09016994374947424102293), "   log(phi^-5) = ", -log(11.09016994374947424102293));
default(realprecision,2600);
\\ control: Zagier D weight-1 small form with exact alpha = zeta(2)/5
b1 = vector(M+1); b1[1]=0.0; b1[2]=1.0; for(n=1,M-1, b1[n+2]=((11*n^2+11*n+3)*b1[n+1]+n^2*b1[n])/(n+1)^2);
g = zeta(2)/5;
r1 = vector(M+1, i, g*a[i]-b1[i]);
default(realprecision,25);
print("Zagier D weight-1 control (alpha = zeta(2)/5 exactly): log|r_n/r_{n-1}| at n=100,200,300,400 : ", log(abs(r1[101]/r1[100])), " ", log(abs(r1[201]/r1[200])), " ", log(abs(r1[301]/r1[300])), " ", log(abs(r1[401]/r1[400])));
\\ and the SQUARE of the weight-1 small form, expanded in the Sym^2 basis
default(realprecision,2600);
q = vector(M+1, i, sum(j=0,i-1, r1[j+1]*r1[i-j]));
default(realprecision,25);
print("[t^n] of (alpha F - G)^2 -- i.e. the genuine Sym^2 small form: log|q_n/q_{n-1}| at n=100,200,300,400 : ", log(abs(q[101]/q[100])), " ", log(abs(q[201]/q[200])), " ", log(abs(q[301]/q[300])), " ", log(abs(q[401]/q[400])));
