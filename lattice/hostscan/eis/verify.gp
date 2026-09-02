default(parisize,1200000000);
cord(f,n) = if(f<=1,1,znorder(Mod(n,f)))
plist(f) = my(G,r); G=znstar(f,1); r=List(); for(n=1,f, if(gcd(n,f)==1, if(type(znconreyconductor(G,n))=="t_INT", listput(r,[n,zncharisodd(G,n),cord(f,n)])))); Vec(r)
PL = vector(60,f,plist(f));
ZS = vector(60,f,znstar(f,1));
indlab(u,lab,N) = if(u==N, lab, znconreyexp(ZS[N], zncharinduce(ZS[u],lab,N)))
dirs(N,k) = my(r,par,aa,bb); r=List(); par=if(k==3,1,0); fordiv(N,u, for(i=1,#PL[u], aa=PL[u][i]; fordiv(N/u,v, for(j=1,#PL[v], bb=PL[v][j]; if((aa[2]+bb[2])%2==par, fordiv(N/(u*v),d, listput(r,[u,aa[1],aa[2],v,bb[1],bb[2],d,lift(Mod(indlab(u,aa[1],N)*indlab(v,bb[1],N),N))]))))))); Vec(r)
chOf(f,lab) = if(f==1, 1, Mod(lab,f))
mkE(k,D) = mfbd(mfeisenstein(k, chOf(D[4],D[5]), chOf(D[1],D[2])), D[7])
ct0e(mf,F) = my(e,pp); e=mfslashexpansion(mf,F,[0,-1;1,0],0,1,&pp); if(pp[1]!=0, return(0)); e[1]
chk(N,k) = my(DD,i,mf,F,a0,b0,p1,p2); DD=dirs(N,k); for(i=1,#DD, mf=if(N==1, mfinit([N,k],3), mfinit([N,k,Mod(DD[i][8],N)],3)); F=mkE(k,DD[i]); a0=mfcoefs(F,0)[1]; b0=ct0e(mf,F); p1=if(a0!=0,1,0); p2=if(b0!=0,1,0); print("CHK ",N," ",k," psi=",DD[i][1],",",DD[i][2]," phi=",DD[i][4],",",DD[i][5]," d=",DD[i][7]," a0nz=",p1," psi1=",if(DD[i][1]==1,1,0)," c0nz=",p2," phi1=",if(DD[i][4]==1,1,0)," scal=",b0*DD[i][7]^k)); 1
GZ=[1,2,3,4,5,6,7,8,9,10,12,13,16,18,25];
for(t=1,#GZ, N=GZ[t]; for(k=3,4, chk(N,k)));
quit;
