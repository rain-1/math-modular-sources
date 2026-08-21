\\ Cooper's three third-order sporadic families s_7, s_10, s_18.
\\ Recurrences (exact, from sporadic_eisenstein_cooper_research_notes.txt lines 346-359):
\\   s_10: (n+1)^3 A_{n+1} = 2(2n+1)(3n^2+3n+1) A_n + 4n(16n^2-1) A_{n-1}
\\   s_7 : (n+1)^3 A_{n+1} = (2n+1)(13n^2+13n+4) A_n + 3n(9n^2-1) A_{n-1}
\\   s_18: (n+1)^3 A_{n+1} = 2(2n+1)(7n^2+7n+3) A_n - 12n(16n^2-1) A_{n-1}
\\ A_0=1, A_1 = actual Cooper IC found empirically by matching archimedean limit
\\ (A_1 = quadratic-constant-term b does NOT work for s_10, s_18; verified against
\\ s_10(n)=sum_k binom(n,k)^4 and against convergence to the claimed limits):
\\   s_10: A_1=2 (matches sum_k binom(n,k)^4: 1,2,18,164,...)
\\   s_7 : A_1=4 (matches F_7 q-series coeffs 1,4,12,16,... in notes)
\\   s_18: A_1=6 (found by scanning A_1=1..8 for convergence to L(2,chi_-3)/2)
default(realprecision,120);
nm=["s_10","s_7","s_18"];
bval=[2,4,6];
N=500;

f10(n,u1,u0)=2*(2*n+1)*(3*n^2+3*n+1)*u1 + 4*n*(16*n^2-1)*u0;
f7(n,u1,u0)=(2*n+1)*(13*n^2+13*n+4)*u1 + 3*n*(9*n^2-1)*u0;
f18(n,u1,u0)=2*(2*n+1)*(7*n^2+7*n+3)*u1 - 12*n*(16*n^2-1)*u0;

rows(fn,b)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=fn(n,A[n+1],A[n])/(n+1)^3;
             B[n+2]=fn(n,B[n+1],B[n])/(n+1)^3);[A,B]};

R10=rows(f10,3); R7=rows(f7,4); R18=rows(f18,6);
R=[R10,R7,R18];

z2=Pi^2/6; Lm3=lfun(lfuncreate(-3),2);
print("=== characteristic roots, c ===");
{my(cr=[[12,-64],[26,-27],[28,192]]);
 for(i=1,3,my(rt=polroots(x^2-cr[i][1]*x+cr[i][2]));
   print(nm[i]," roots(asymp x^2-",cr[i][1],"x+",cr[i][2],")= ",rt~," product=",cr[i][2]));}

print("\n=== limit b_n/a_n, lindep vs [1,zeta(2),L(2,chi-3)] ===");
{for(i=1,3,my(A=R[i][1],B=R[i][2]);
  my(l=B[N+1]/A[N+1]);
  print(nm[i]," lim(n=",N,")=",l);
  print("   lindep[1,z2,L2chi-3]=",lindep([l,1,z2,Lm3])~));}

print("\n=== p-adic slopes p=2,3,5,7, max v_p(A_n) for n<=",N," ===");
{for(i=1,3,my(A=R[i][1],B=R[i][2]);
 my(sl=vector(4,j,my(p=[2,3,5,7][j]); valuation(B[N+1]/A[N+1]-B[N]/A[N],p)*1./N));
 my(mv=vector(4,j,my(p=[2,3,5,7][j]); vecmax(vector(N,n,valuation(A[n+1],p)))));
 print(nm[i]," slopes=",sl,"  maxv_p(A_n)=",mv));}

print("\n=== denominator check: is (n+1)*A_n an integer, and lcm(1..n)^2*B_n integral? ===");
{for(i=1,3,my(A=R[i][1],B=R[i][2]);
 my(dA=vector(30,n,denominator(A[n+1]))); \\ n=0..29
 my(dB=vector(30,n,denominator(lcm(vector(n,k,k))^2*B[n+1])));
 print(nm[i]," denominators of A_n, n=0..29: ",dA);
 print(nm[i]," denom of lcm(1..n)^2*B_n, n=1..29 (want all 1): ",vecextract(dB,"2..30")));}
\q
