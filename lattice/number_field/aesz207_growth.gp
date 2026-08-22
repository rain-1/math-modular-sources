default(parisizemax, 4000000000);
default(realprecision, 2600);
Q0(t)=t^4
Q1(t)=-2^4*(1072*t^4 - 17824*t^3 - 10888*t^2 - 1976*t - 145)
Q2(t)=-2^17*(51088*t^4 + 116368*t^3 - 45264*t^2 - 14228*t - 1397)
Q3(t)=2^28*13*(73104*t^4 + 1536*t^3 - 488*t^2 + 384*t + 97)
Q4(t)=-2^44*13^2*(2*t+1)^4
fwd(v0,N) = my(u=vector(N+1)); for(i=1,4,u[i]=v0[i]*1.0); for(n=4,N, u[n+1] = -( Q1(n-1)*u[n] + Q2(n-2)*u[n-1] + Q3(n-3)*u[n-2] + Q4(n-4)*u[n-3] )/Q0(n)); u
N=520;
print("--- forward growth from the 4 coordinate initial vectors (u_0..u_3) ---");
for(i=1,4, my(v=vector(4,j,if(j==i,1,0)), u=fwd(v,N)); print("  e_",i," : log|u_n|/n at n=",N," = ", log(abs(u[N+1]))/N, "   log|u_n/u_{n-1}| = ", log(abs(u[N+1]/u[N]))));
print("  target log 89531.389... = ", log(89531.3892067201467820569392329979580336015275720649903262571));
print("  log 53248 = ", log(53248.0), "   log 187.389... = ", log(187.389206720146782056939232997958033601527572064990326257133));
print("");
print("--- backward recurrence: dominant backward solution = smallest forward growth ---");
NB=520;
bkw(v0) = my(u=vector(NB+1)); for(i=0,3, u[NB+1-i]=v0[4-i]*1.0); forstep(n=NB,4,-1, u[n-3] = -( Q0(n)*u[n+1] + Q1(n-1)*u[n] + Q2(n-2)*u[n-1] + Q3(n-3)*u[n-2] )/Q4(n-4)); u
w1 = bkw([1,0,0,0]); w2 = bkw([0,1,0,0]); w3 = bkw([0,0,1,0]); w4 = bkw([0,0,0,1]);
v1 = vector(4,j,w1[j]/w1[1]); v2 = vector(4,j,w2[j]/w2[1]); v3 = vector(4,j,w3[j]/w3[1]);
print("  normalised recovered initial vectors (u_0=1) from 3 different backward seeds:");
default(realprecision,40);
print("   seed e1 : ", v1);
print("   seed e2 : ", v2);
print("   seed e3 : ", v3);
default(realprecision,2600);
print("  agreement seed1 vs seed2, digits: ", vector(4,j,if(v1[j]==v2[j],"exact",round(-log(abs(v1[j]-v2[j])+1e-2500)/log(10)))));
print("  agreement seed1 vs seed3, digits: ", vector(4,j,if(v1[j]==v3[j],"exact",round(-log(abs(v1[j]-v3[j])+1e-2500)/log(10)))));
uu = fwd(v1, 400);
print("  forward run from the recovered vector: log|u_n|/n at n=100,200,300,400 = ");
print("     ", log(abs(uu[101]))/100, " / ", log(abs(uu[201]))/200, " / ", log(abs(uu[301]))/300, " / ", log(abs(uu[401]))/400);
print("  consecutive ratio log|u_400/u_399| = ", log(abs(uu[401]/uu[400])));
print("  log 187.3892067... = ", log(187.389206720146782056939232997958033601527572064990326257133));
print("");
print("--- is the subdominant initial vector defined over Q(sqrt17) (or Q)? ---");
default(realprecision, 300);
s17 = sqrt(17);
for(j=2,4, print("  u_",j-1,"/u_0 = ", v1[j]));
for(j=2,4, print("   lindep([1, sqrt17, u_",j-1,"/u_0])   (200 digits) = ", lindep(vector(3,i,[1,s17,v1[j]][i]), 200)));
for(j=2,4, print("   lindep([1, u_",j-1,"/u_0])           (200 digits) = ", lindep([1,v1[j]], 200)));
print("");
print("--- the 'kill the dominant' functional: ratios u^{(i)}_N/u^{(1)}_N ---");
default(realprecision, 2600);
f1 = fwd([1,0,0,0],N); f2 = fwd([0,1,0,0],N); f3 = fwd([0,0,1,0],N); f4 = fwd([0,0,0,1],N);
rr = [f2[N+1]/f1[N+1], f3[N+1]/f1[N+1], f4[N+1]/f1[N+1]];
default(realprecision, 300);
for(j=1,3, print("   r_",j," = ", rr[j]));
for(j=1,3, print("   lindep([1, sqrt17, r_",j,"]) = ", lindep([1,s17,rr[j]], 200)));
