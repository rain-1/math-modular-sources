default(parisize,900000000);
ct0(mf,f) = my(e,pp); e=mfslashexpansion(mf,f,[0,-1;1,0],0,1,&pp); if(pp[1]!=0, return([0,pp[1],pp[3]])); [e[1],pp[1],pp[3]]
doann(N,k,m,mf) = my(B,d,r1,r2,M,K,flag); B=mfbasis(mf); d=#B; r1=vector(d,i,mfcoefs(B[i],0)[1]); r2=vector(d,i,ct0(mf,B[i])[1]); flag=vector(d,i,ct0(mf,B[i])[2]); M=matconcat([Col(r1),Col(r2)]~); K=matker(M); print("ANN ",N," ",k," ",m," ",d," ",#K," ",vecmax(flag)," ",d-matrank(M))
GZ=[1,2,3,4,5,6,7,8,9,10,12,13,16,18,25];
for(t=1,#GZ, N=GZ[t]; for(k=3,4, for(m=1,N, if(gcd(m,N)==1, if(zncharisodd(znstar(N,1),m)==(k%2), doann(N,k,m, if(N==1, mfinit([N,k],3), mfinit([N,k,Mod(m,N)],3))))))));
quit;
