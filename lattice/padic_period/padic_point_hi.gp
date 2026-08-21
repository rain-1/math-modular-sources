default(parisize,6000000000);
N=1800; M=1600; 
rows(a,b,c,n)={my(A=vector(n+1),B=vector(n+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(m=1,n-1,A[m+2]=((a*m^2+a*m+b)*A[m+1]-c*m^2*A[m])/(m+1)^2;
            B[m+2]=((a*m^2+a*m+b)*B[m+1]-c*m^2*B[m])/(m+1)^2);[A,B]};
EV(c,t,K)={my(s=O(3^M),T=1+O(3^M)); for(k=0,K,s+=c[k+1]*T;T*=t); s};
gettime(); RC=rows(10,3,9,N); print("rows time ",gettime());
xiQ = RC[2][N+1]/RC[1][N+1]; xi = xiQ + O(3^M);
sol(p,K)={my(x0=p[1],u=p[2],v=p[3],ub=p[4],vb=p[5],G,Gb,SA,SB,c0,c1,SAf,SBf);
 G=(1+3*v)/(1-3*u); Gb=(1+3*vb)/(1-3*ub);
 c1=(G*EV(RC[1],u,K)+G^2*EV(RC[1],v,K)/4)-(Gb*EV(RC[1],ub,K)+Gb^2*EV(RC[1],vb,K)/4);
 c0=(G*EV(RC[2],u,K)+G^2*EV(RC[2],v,K)/4)-(Gb*EV(RC[2],ub,K)+Gb^2*EV(RC[2],vb,K)/4);
 print("   K=",K," v3(c0)=",valuation(c0)," v3(c1)=",valuation(c1)," v3(c0/c1 - xi)=",valuation(c0/c1-xi));
};
p1=[-1,3/2,1/4,-3/5,1/25]; p2=[-3,15/7,9/49,-5/4,1/16]; p3=[-6,14/5,4/25,-21/11,9/121];
print("point x0=-1"); sol(p1,300); sol(p1,600); sol(p1,900);
print("point x0=-3"); sol(p2,900);
print("point x0=-6"); sol(p3,900);
\\ direct residual with full xi
res(p,K)={my(x0=p[1],u=p[2],v=p[3],ub=p[4],vb=p[5],G,Gb,r,Rf,H,Hb);
 r=vector(K+1,k,RC[2][k]-xi*RC[1][k]);
 G=(1+3*v)/(1-3*u); Gb=(1+3*vb)/(1-3*ub);
 H=G*EV(r,u,K)+G^2*EV(r,v,K)/4; Hb=Gb*EV(r,ub,K)+Gb^2*EV(r,vb,K)/4;
 print("   K=",K," v3(H)=",valuation(H)," v3(H-Hb)=",valuation(H-Hb));
};
print("residuals"); res(p1,900); res(p2,900); res(p3,900);
quit
