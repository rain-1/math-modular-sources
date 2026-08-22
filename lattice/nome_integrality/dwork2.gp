default(parisizemax,6000000000);
M=400;
A=vector(M+2); A[1]=1; A[2]=5;
{for(n=1,M, A[n+2] = ((2*n+1)*(17*n^2+17*n+5)*A[n+1] - n^3*A[n])/(n+1)^3);}
B=vector(M+2); B[1]=0;
{for(n=1,M+1, my(m=n-1,bm1,bm2,am,am1,am2); bm1=B[n]; bm2=if(n>=2,B[n-1],0); am=A[n+1]; am1=A[n]; am2=if(n>=2,A[n-1],0); B[n+1] = ( (2*n-1)*(17*m^2+17*m+5)*bm1 - m^3*bm2 - 3*n^2*am + (102*m^2+102*m+27)*am1 - 3*m^2*am2 )/n^3 );}
print("B_1..B_6 = ", vector(6,i,B[i+1]));
print("B_59/A_59 = ", B[60]*1.0/A[60], "   2*zeta(3) = ", 2*zeta(3));
Fs = sum(n=0,M, A[n+1]*t^n) + O(t^(M+1));
Gs = sum(n=0,M, B[n+1]*t^n) + O(t^(M+1));
sig = sqrt(1-34*t+t^2+O(t^(M+1)));
uu = Gs/Fs;
print("check theta(G/F) == 1/(F sigma) - 1 : ", t*deriv(uu) - (1/(Fs*sig) - 1));
quit;
