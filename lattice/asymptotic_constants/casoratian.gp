default(parisize,2000000000);
default(realprecision,400);
\\ rows: [name, type(2 or 3), a, b, c, d, a1, N, lambda1, lambda2, kappa]
rows=[["Apery",3,17,5,1,0,5,6],["Domb",3,10,4,64,0,4,12],["eps",3,12,4,16,0,4,8],["zeta",3,9,3,-27,0,3,9],["s7",3,13,4,-27,3,4,7],["s10",3,6,2,-64,4,2,10],["s18",3,14,6,192,-12,6,18],["FranelA",2,7,2,-8,0,2,0]];
P3(n,a,b)=(2*n+1)*(a*n^2+a*n+b);
Q3(n,c,d)=n*(c*n^2+d);
P2(n,a,b)=a*n^2+a*n+b;
Q2(n,c)=c*n^2;
{foreach(rows,R, my(nm=R[1],ty=R[2],a=R[3],b=R[4],c=R[5],d=R[6],a1=R[7],N=R[8]);
  my(l1,l2); if(ty==3, l1=a+sqrt(a^2-c); l2=a-sqrt(a^2-c), l1=(a+sqrt(a^2-4*c))/2; l2=(a-sqrt(a^2-4*c))/2);
  my(s=ty/2);
  \\ forward: a_n exact to M
  my(M=1500, A=vector(M+2)); A[1]=1; A[2]=a1;
  for(n=1,M, if(ty==3, A[n+2]=(P3(n,a,b)*A[n+1]-Q3(n,c,d)*A[n])/(n+1)^3, A[n+2]=(P2(n,a,b)*A[n+1]-Q2(n,c)*A[n])/(n+1)^2));
  \\ backward (Miller) recessive solution from N0
  my(N0=6000, W=vector(N0+2)); W[N0+2]=0.; W[N0+1]=1.;
  forstep(n=N0,1,-1, if(ty==3, W[n]=(P3(n,a,b)*W[n+1]-(n+1)^3*W[n+2])/Q3(n,c,d), W[n]=(P2(n,a,b)*W[n+1]-(n+1)^2*W[n+2])/Q2(n,c))); 
  \\ W[n+1] ~ w_n ; normalise a_1 w_0 - a_0 w_1 = 1
  my(nrm=a1*W[1]-W[2]); W=W/nrm;
  my(xi=W[1]);
  \\ Richardson for K+ and K-
  my(rich=(f,m,Nn,st)->my(Mt=matrix(m,m),v=vector(m)); for(i=1,m,my(n=Nn-st*(i-1)); for(j=1,m,Mt[i,j]=1.0/n^(j-1)); v[i]=f(n)); (matsolve(Mt,v~))[1]);
  my(Kp=rich(n->A[n+1]*1.0/l1^n*n^s,20,M,50));
  my(Km=rich(n->W[n+1]/l2^n*n^s,20,N0-1000,100));
  my(kap); if(ty==3&&d!=0, my(r=d/c); if(r>0, kap=sinh(Pi*sqrt(r))/(Pi*sqrt(r)), kap=sin(Pi*sqrt(-r))/(Pi*sqrt(-r))), kap=1);
  my(Kp_th); if(ty==3, Kp_th=sqrt(N)/(2*Pi^(3/2))*sqrt(l1/(l1-l2)), Kp_th=0);
  print(nm,": xi=",xi*1.0+0.," K+=",Kp," K+_thm=",Kp_th,"  K-=",Km,"  K+K-(l1-l2)=",Kp*Km*(l1-l2),"  kappa=",kap,"  K-_pred=",if(ty==3,kap/((l1-l2)*Kp_th),0));
);}
print("Apery classical (b1=6): K- = ", 6*2^(-1/4)/6*Pi^(3/2)*(sqrt(2)-1)^2, " vs 6*K-(b1=1)");
print("predictions Cooper: s7 sqrt(pi)/14=",sqrt(Pi)/14,"  s10 sqrt(pi)/10=",sqrt(Pi)/10,"  s18 sqrt(pi)/6=",sqrt(Pi)/6,"  Domb pi^(3/2)/24=",Pi^(3/2)/24);
quit;
