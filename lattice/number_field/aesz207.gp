default(parisizemax, 4000000000);
default(realprecision, 60);
Q0(t)=t^4
Q1(t)=-2^4*(1072*t^4 - 17824*t^3 - 10888*t^2 - 1976*t - 145)
Q2(t)=-2^17*(51088*t^4 + 116368*t^3 - 45264*t^2 - 14228*t - 1397)
Q3(t)=2^28*13*(73104*t^4 + 1536*t^3 - 488*t^2 + 384*t + 97)
Q4(t)=-2^44*13^2*(2*t+1)^4
print("=== 14. AESZ 207");
\\ --- the holomorphic solution, recurrence imposed for all n>=1 -------------
N=200;
u = vector(N+1); u[1]=1;
for(n=1,N, u[n+1] = -( if(n>=1,Q1(n-1)*u[n],0) + if(n>=2,Q2(n-2)*u[n-1],0) + if(n>=3,Q3(n-3)*u[n-2],0) + if(n>=4,Q4(n-4)*u[n-3],0) )/Q0(n));
print("  A_0..A_4 : ", vector(5,i,u[i]));
print("  A_1 = -2320 ? ", u[2]==-2320, "    A_2 = 57601296 ? ", u[3]==57601296);
print("  A_n integral for n<=",N," ? ", if(sum(n=0,N,denominator(u[n+1])!=1)==0,"yes","NO"));
\\ --- characteristic polynomial (leading coefficients) ---------------------
c0=1; c1=-2^4*1072; c2=-2^17*51088; c3=2^28*13*73104; c4=-2^44*13^2*16;
cp = c0*x^4+c1*x^3+c2*x^2+c3*x+c4;
print("  leading coeffs c0..c4 = ", [c0,c1,c2,c3,c4]);
print("  charpoly = ", cp);
tst = (x-53248)^2*(x^2+89344*x-2^24);
print("  claimed factorisation (x-53248)^2 (x^2+89344x-2^24) = ", tst);
print("  EXACT MATCH ? ", cp == tst);
print("  factor over Q : ", factor(cp));
print("  53248 = 2^12*13 ? ", 53248 == 2^12*13);
print("  roots : ", polroots(cp));
r1 = -2^7*(349+85*sqrt(17)); r2 = -2^7*(349-85*sqrt(17));
print("  -2^7(349+85 sqrt17) = ", r1, "    -2^7(349-85 sqrt17) = ", r2);
print("  quadratic factor disc = ", poldisc(x^2+89344*x-2^24), " = 2^?*17: ", factor(poldisc(x^2+89344*x-2^24)));
print("  N_{Q(sqrt17)/Q} of 2^7(349+85 sqrt17) = ", 2^14*(349^2-85^2*17), "  |N| = ", abs(2^14*(349^2-85^2*17)), " = 2^24 ? ", abs(2^14*(349^2-85^2*17))==2^24);
print("  |N|^(1/2) = ", sqrt(abs(2^14*(349^2-85^2*17))), " = 4096 = 2^12");
\\ --- margins --------------------------------------------------------------
kk=4;
print("  (i)  row over Q, one archimedean place; lambda_1 = ", abs(r1), ", lambda_2 = 53248 (double)");
print("       rho = 1/53248 = ", 1.0/53248, "  log rho = ", log(1.0/53248), "  margin (k=4) = ", log(1.0/53248)-kk);
print("  (ii) counterfactual Q(sqrt17) 2-place geometric mean = ", sqrt(abs(r1*r2)), " = 4096");
print("       log(4096) = ", log(4096.0), "  margin = log4096 - 4 = ", log(4096.0)-kk);
print("       12*log2 = ", 12*log(2.0));
\\ --- companion and sharp k -------------------------------------------------
NK=150;
LCT = vector(NK+2); LCT[1]=1; for(n=1,NK+1, LCT[n+1]=lcm(LCT[n],n));
b = vector(N+1); b[1]=0; b[2]=1; b[3]=0; b[4]=0;
for(n=4,N, b[n+1] = -( Q1(n-1)*b[n] + Q2(n-2)*b[n-1] + Q3(n-3)*b[n-2] + Q4(n-4)*b[n-3] )/Q0(n));
badat(v,k,M) = my(bad=-1); for(n=1,M, if(denominator(LCT[n+1]^k*v[n+1])!=1, bad=n; break)); bad
sharpk(v,M) = my(kk2=-1); for(k=0,10, if(badat(v,k,M)<0, kk2=k; break)); kk2
k14 = sharpk(b,NK);
print("  companion (0,1,0,0), recurrence imposed n>=4: sharp k (n<=",NK,") = ", k14, "  first fail at k-1: n=", badat(b,k14-1,NK));
