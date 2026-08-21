\\ Cross-row alignment: Cooper s_10,s_7,s_18 vs Zagier second-order rows sharing a period.
default(realprecision,250);
f10(n,u1,u0)=2*(2*n+1)*(3*n^2+3*n+1)*u1 + 4*n*(16*n^2-1)*u0;
f7(n,u1,u0)=(2*n+1)*(13*n^2+13*n+4)*u1 + 3*n*(9*n^2-1)*u0;
f18(n,u1,u0)=2*(2*n+1)*(7*n^2+7*n+3)*u1 - 12*n*(16*n^2-1)*u0;
rowsC(fn,b,N)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=fn(n,A[n+1],A[n])/(n+1)^3;
             B[n+2]=fn(n,B[n+1],B[n])/(n+1)^3);[A,B]};
rowsZ(a,b,c,N)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;
            B[n+2]=((a*n^2+a*n+b)*B[n+1]-c*n^2*B[n])/(n+1)^2);[A,B]};
N=300;
R10=rowsC(f10,2,N); R7=rowsC(f7,4,N); R18=rowsC(f18,6,N);
RA=rowsZ(7,2,-8,N); RD=rowsZ(11,3,-1,N); RC=rowsZ(10,3,9,N);
RB=rowsZ(9,3,27,N); RF=rowsZ(17,6,72,N);

crossdet(r,rp,X,Y,ns)={vector(#ns,k,my(n=ns[k]);valuation(r*X[1][n+1]*Y[2][n+1]-rp*Y[1][n+1]*X[2][n+1],p_used))};

pairtest(name,r,rp,X,Y,p)={my(ns=[50,100,150,200,250,300]);
 print(name," p=",p,": v_p at n=",ns," = ",
  vector(#ns,k,my(n=ns[k]);valuation(r*X[1][n+1]*Y[2][n+1]-rp*Y[1][n+1]*X[2][n+1],p)))};

print("=== s_10 (lim z2/5) vs Zagier D (lim z2/5): r=r'=1 ===");
{for(pp=1,4,pairtest("s10,D",1,1,R10,RD,[2,3,5,7][pp]));}

print("\n=== s_7 (lim z2/7) vs Zagier A (lim z2/4): r=7,r'=4 so both ->z2 ===");
{for(pp=1,4,pairtest("s7,A",7,4,R7,RA,[2,3,5,7][pp]));}

print("\n=== s_18 (lim L/2) vs Zagier C (lim L/2): r=r'=1 ===");
{for(pp=1,4,pairtest("s18,C",1,1,R18,RC,[2,3,5,7][pp]));}

print("\n=== s_18 (lim L/2) vs Zagier B (lim L/2, no arch limit): r=r'=1 ===");
{for(pp=1,4,pairtest("s18,B",1,1,R18,RB,[2,3,5,7][pp]));}

print("\n=== s_18 (lim L/2) vs Zagier F (lim 5L/8): r=5,r'=4 so both ->L ===");
{for(pp=1,4,pairtest("s18,F",5,4,R18,RF,[2,3,5,7][pp]));}
\q
