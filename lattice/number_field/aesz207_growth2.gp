default(parisizemax, 4000000000);
default(realprecision, 2600);
Q0(t)=t^4
Q1(t)=-2^4*(1072*t^4 - 17824*t^3 - 10888*t^2 - 1976*t - 145)
Q2(t)=-2^17*(51088*t^4 + 116368*t^3 - 45264*t^2 - 14228*t - 1397)
Q3(t)=2^28*13*(73104*t^4 + 1536*t^3 - 488*t^2 + 384*t + 97)
Q4(t)=-2^44*13^2*(2*t+1)^4
fwd(v0,N) = my(u=vector(N+1)); for(i=1,4,u[i]=v0[i]*1.0); for(n=4,N, u[n+1] = -( Q1(n-1)*u[n] + Q2(n-2)*u[n-1] + Q3(n-3)*u[n-2] + Q4(n-4)*u[n-3] )/Q0(n)); u
N=520;
f = [fwd([1,0,0,0],N), fwd([0,1,0,0],N), fwd([0,0,1,0],N), fwd([0,0,0,1],N)];
default(realprecision,25);
print("consecutive-ratio growth log|u_n/u_{n-1}| for the 4 coordinate seeds:");
for(i=1,4, print("  e_",i,":  n=200: ",log(abs(f[i][201]/f[i][200])), "   n=400: ", log(abs(f[i][401]/f[i][400])), "   n=520: ", log(abs(f[i][521]/f[i][520]))));
print("  log(89531.3892067201468) = ", log(89531.3892067201467820569392329979580336015275720649903));
default(realprecision,2600);
NB=520;
bkw(v0) = my(u=vector(NB+1)); for(i=0,3, u[NB+1-i]=v0[4-i]*1.0); forstep(n=NB,4,-1, u[n-3] = -( Q0(n)*u[n+1] + Q1(n-1)*u[n] + Q2(n-2)*u[n-1] + Q3(n-3)*u[n-2] )/Q4(n-4)); u
w1 = bkw([1,0,0,0]);
v1 = vector(4,j,w1[j]/w1[1]);
uu = fwd(v1, 460);
default(realprecision,25);
print("subdominant solution recovered by backward iteration, initial vector (u_0=1):");
print("  ", v1);
print("  log|u_n|/n            : n=100 ",log(abs(uu[101]))/100,"  n=200 ",log(abs(uu[201]))/200,"  n=300 ",log(abs(uu[301]))/300,"  n=400 ",log(abs(uu[401]))/400);
print("  log|u_n/u_{n-1}|      : n=100 ",log(abs(uu[101]/uu[100])),"  n=200 ",log(abs(uu[201]/uu[200])),"  n=300 ",log(abs(uu[301]/uu[300])),"  n=400 ",log(abs(uu[401]/uu[400])),"  n=460 ",log(abs(uu[461]/uu[460])));
print("  log(187.3892067201468) = ", log(187.389206720146782056939232997958033601527572064990326257133));
print("  log(53248) = ", log(53248.0));
