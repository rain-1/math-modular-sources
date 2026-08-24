default(parisize, 4000000000);
default(realprecision, 60);
{P6 = [nn^2+2*nn+1, 9*nn^2+nn, 20*nn^2-40*nn+20, -20*nn^2+20*nn+4, -96*nn^2+352*nn-336, -64*nn^2+320*nn-400];}
NMX = 40000;
{go(u1) = my(r=5, U=vector(NMX+1));
  if(u1==0, U[1]=1.0, U[2]=1.0);
  if(u1==0, U[2]=0.0);
  /* placement 6 modular: a_0=1,a_1=0,a_2=0,a_3=4,a_4=-12,a_5=48 */
  if(u1==0, U[1]=1.0;U[2]=0.0;U[3]=0.0;U[4]=4.0;U[5]=-12.0;U[6]=48.0);
  my(st = if(u1==0, 5, 1));
  for(n=st, NMX-1, my(s=0.0);
    for(j=1,r, my(idx=n+1-j); if(idx>=0, s += subst(P6[j+1],nn,n)*U[idx+1]));
    U[n+2] = -s/subst(P6[1],nn,n));
  U;}
A = go(0); W = go(1);
print("a_n check: ", vector(10,i,round(A[i])));
print("w_n check: ", vector(6,i,W[i]));
NS = [500,1000,2000,4000,8000,16000,32000,40000];
{R = vector(#NS, i, W[NS[i]+1]/A[NS[i]+1]);}
for(i=1,#NS, print("  n=",NS[i],"  L=log n=",log(NS[i]*1.0),"  xi(n) = ", R[i]));
/* fit xi + c1/L + c2/L^2 + ... using the last k points */
{for(k=2,6,
   my(idx = vector(k,i,#NS-k+i));
   my(MM = matrix(k,k,i,j, my(Lg=log(NS[idx[i]]*1.0)); if(j==1,1.0,Lg^(-(j-1)))));
   my(rr = vectorv(k,i,R[idx[i]]));
   my(sol = matsolve(MM,rr));
   print("  fit with ",k," terms (n from ",NS[idx[1]],"): xi_inf = ", sol[1]));}
print("G/9 = ", Catalan/9, "   (2/9)zeta(2) = ", 2*zeta(2)/9, "   1 = 1");
quit;
