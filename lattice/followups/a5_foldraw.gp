default(parisize, 6000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
default(realprecision, 40);
O184=0; for(i=1,#OPS, if(OPS[i][1]=="184", O184=OPS[i]));
rowR3(a,b,c,N)={my(A=vector(N+1),B=vector(N+1)); A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1, A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
              B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3); [A,B];}
N=2500;
pe=rowR3(11,5,125,N); Ae=vector(N+1,i,1.0*pe[1][i]); Be=vector(N+1,i,1.0*pe[2][i]);
pr=aperyPair(O184[4],N); A4=vector(N+1,i,1.0*pr[1][i]); B4=vector(N+1,i,1.0*pr[2][i]);
L53=lfun(5,3); L52=lfun(5,2);
tgtE = L53/2 + I*Pi*L52/10;  tgt4 = L53/4 + I*Pi*L52/20;
xin(A,B,rho,n) = (B[n+2]-rho*B[n+1])/(A[n+2]-rho*A[n+1]);
print("target eta = ", tgtE);
print("target 184 = ", tgt4);
print("\nraw xi_n (eta, rho=11-2I):");
forstep(n=2000,2400,100, my(x=xin(Ae,Be,11-2*I,n)); print("  n=",n," ",x,"  |x-tgtE|=",abs(x-tgtE)," |x-conj(tgtE)|=",abs(x-conj(tgtE))));
print("\nraw xi_n (eta, rho=11+2I):");
forstep(n=2000,2400,100, my(x=xin(Ae,Be,11+2*I,n)); print("  n=",n," ",x,"  |x-tgtE|=",abs(x-tgtE)," |x-conj(tgtE)|=",abs(x-conj(tgtE))));
print("\nwindow avg n=1000..2400 (eta, rho=11-2I):");
{my(s=0.); for(n=1000,2400,s+=xin(Ae,Be,11-2*I,n)); s/=1401; print("  ",s,"   diff to tgtE = ",abs(s-tgtE),"  diff to conj = ",abs(s-conj(tgtE)));}
print("\nwindow avg n=1000..2400 (184, rho=44-8I):");
{my(s=0.); for(n=1000,2400,s+=xin(A4,B4,44-8*I,n)); s/=1401; print("  ",s,"   diff to tgt4 = ",abs(s-tgt4),"  diff to conj = ",abs(s-conj(tgt4)));}
print("\nratio xi_n(184)/xi_n(eta) at n=2400: ", xin(A4,B4,44-8*I,2400)/xin(Ae,Be,11-2*I,2400));
quit
