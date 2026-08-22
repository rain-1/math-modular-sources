default(parisizemax, 6000000000);
default(realprecision, 60);
\\ K = Q(sqrt5), w=(1+sqrt5)/2, w^2=w+1
ww = Mod(y, y^2-y-1);
s5 = 2*ww-1;
co(e) = my(p=lift(e)); [polcoeff(p,0), polcoeff(p,1)];
isint(e) = my(v=co(e)); (denominator(v[1])==1)&&(denominator(v[2])==1);
tos5(e) = my(v=co(e)); [v[1]+v[2]/2, v[2]/2];  \\ a + b*sqrt5
AA = 124+55*s5; BB = 34+15*s5; CC = 1;
Pf(n) = (2*n+1)*(AA*n^2+AA*n+BB);
Qf(n) = CC*n^3;
N = 220;
u = vector(N+2); u[1]=1; u[2]=Pf(0)*u[1];
for(n=1,N, u[n+2] = (Pf(n)*u[n+1]-Qf(n)*u[n])/(n+1)^3);
print("u_1 = ", tos5(u[2]), "  target [34,15]  -> ", u[2]==34+15*s5);
print("u_2 = ", tos5(u[3]), "  target [7111,3180] -> ", u[3]==7111+3180*s5);
print("u_3 = ", tos5(u[4]), "  target [2040334,912465] -> ", u[4]==2040334+912465*s5);
bad=-1; for(n=0,200, if(!isint(u[n+1]), bad=n; break));
print("first non-integral u_n, n<=200: ", if(bad<0,"NONE (all in O_K)",bad));
bad=-1; for(n=0,N, if(!isint(u[n+1]), bad=n; break));
print("first non-integral u_n, n<=",N,": ", if(bad<0,"NONE",bad));
print("first ten u_n as [a,b] meaning a+b*sqrt5:");
for(j=0,9, print("  u_",j," = ", tos5(u[j+1])));
\\ ---- companion + sharp k ----
LCT=vector(200); LCT[1]=1; for(n=1,199, LCT[n+1]=lcm(LCT[n],n));
b = vector(140); b[1]=0; b[2]=1;
for(n=1,138, b[n+2] = (Pf(n)*b[n+1]-Qf(n)*b[n])/(n+1)^3);
okk(k,M) = my(f=-1); for(n=1,M, if(!isint(LCT[n+1]^k*b[n+1]), f=n; break)); f;
for(k=0,5, my(f=okk(k,120)); print("  k=",k," : ", if(f<0,"ALL INTEGRAL to n=120", Str("fails first at n=",f))));
\\ ---- roots and margins ----
default(realprecision, 60);
ev(e,sg) = my(v=co(e)); v[1]+v[2]*(1+sg*sqrt(5))/2;
mar(A,C,k) = {
  my(rts=vector(2), lam1=vector(2), lam2=vector(2));
  for(i=1,2, my(sg=if(i==1,1,-1), a=ev(A,sg), c=ev(C,sg), d=a^2-c, r1,r2);
     if(d>=0, r1=abs(a+sqrt(d)); r2=abs(a-sqrt(d)), r1=sqrt(abs(c)); r2=r1);
     if(r1<r2, my(t=r1); r1=r2; r2=t);
     lam1[i]=r1; lam2[i]=r2);
  print("  place v1 (sqrt5->+): lambda1=",lam1[1]," lambda2=",lam2[1]);
  print("  place v2 (sqrt5->-): lambda1=",lam1[2]," lambda2=",lam2[2]);
  for(v0=1,2, my(vo=3-v0, rho0=1/lam2[v0], rhoo=1/lam1[vo], gm=sqrt(rho0*rhoo));
     print("  v0=v",v0,": rho_v0=",rho0," rho_other=",rhoo," geom.mean=",gm," ModeI margin=",log(gm)-k));
  my(nn=lam2[1]*lam2[2]);
  print("  N(lambda2)=",nn,"  sqrt=",sqrt(nn),"  ModeII margin=",-log(nn)/2-k);
};
mar(AA,CC,3);
quit;
