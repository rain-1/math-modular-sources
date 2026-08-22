N=40;
A=vector(N+2); A[1]=1; A[2]=5;
{for(n=1,N, A[n+2] = ((2*n+1)*(17*n^2+17*n+5)*A[n+1] - n^3*A[n])/(n+1)^3);}
Fs = sum(n=0,N, A[n+1]*t^n) + O(t^(N+1));
R = 1/(Fs*sqrt(1-34*t+t^2+O(t^(N+1))));
r = vector(N+1,i,polcoeff(R,i-1));
c = vector(N); {for(m=1,N, c[m] = (r[m+1] - sum(d=1,m-1, if(m%d==0, d*c[d], 0)))/m);}
print("c_1..c_10 = ", vector(10,i,c[i]));
b=0; {for(m=1,N, if(denominator(c[m])!=1,b++));} print("c_n integral to n=",N,"? viol=",b);
print("A_0..A_6 = ", vector(7,i,A[i]));
quit;
