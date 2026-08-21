default(parisizemax,2000000000); default(realprecision,120);
N=400; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=10;B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=((136*n^2+68*n+10)*A[n+1]-4*(2*n-1)^2*A[n])/(n+1)^2;B[n+2]=((136*n^2+68*n+10)*B[n+1]-4*(2*n-1)^2*B[n])/(n+1)^2);
lim=B[N+1]*1./A[N+1];
for(D=1,30, if(issquarefree(D), print("pi^2/sqrt(",D,"): ",lindep([lim,Pi^2/sqrt(D),1])~)));
for(D=1,30, if(issquarefree(D), print("pi^2*sqrt(",D,")/  lindep: ",lindep([lim,Pi^2*sqrt(D)])~)));
\\ cusp forms: weight 3 newforms of level N with odd quadratic character
{for(Nl=6,48, forstep(D=-1,-30,-1, if(issquarefree(abs(D))||abs(D)%4==0, my(chi); if(D==-4,chi=-4,D==-8,chi=-8,D==-3,chi=-3,D==-24,chi=-24,D==-20,chi=-20,D==-7,chi=-7,1,next);
 if(Nl%abs(chi)!=0,next); my(mf=mfinit([Nl,3,chi],0)); if(mfdim(mf)==0,next);
 my(Lb=mfeigenbasis(mf)); for(i=1,#Lb, my(F=Lb[i]); if(poldegree(mfparams(F)[4])>1,next); my(L=lfunmf(mf,F)); my(val=lfun(L,2)); my(r=lindep([lim,val,1]));
 if(abs(r[3])<2&&abs(r[1])<200&&abs(r[2])<200, print("HIT N=",Nl," chi=",chi," form ",i," : ",r~)));)))}
\q
