default(parisize,1200000000);
default(realprecision,60);
cord(f,n) = if(f<=1,1,znorder(Mod(n,f)))
plist(f) = my(G,r); G=znstar(f,1); r=List(); for(n=1,f, if(gcd(n,f)==1, if(type(znconreyconductor(G,n))=="t_INT", listput(r,[n,zncharisodd(G,n),cord(f,n)])))); Vec(r)
PL = vector(60,f,plist(f));
ZS = vector(60,f,znstar(f,1));
indlab(u,lab,N) = if(u==N, lab, znconreyexp(ZS[N], zncharinduce(ZS[u],lab,N)))
dirs(N,k) = my(r,par,aa,bb); r=List(); par=if(k==3,1,0); fordiv(N,u, for(i=1,#PL[u], aa=PL[u][i]; fordiv(N/u,v, for(j=1,#PL[v], bb=PL[v][j]; if((aa[2]+bb[2])%2==par, fordiv(N/(u*v),d, listput(r,[u,aa[1],aa[2],v,bb[1],bb[2],d,lift(Mod(indlab(u,aa[1],N)*indlab(v,bb[1],N),N))]))))))); Vec(r)
chOf(f,lab) = if(f==1, 1, Mod(lab,f))
mkE(k,D) = mfbd(mfeisenstein(k, chOf(D[4],D[5]), chOf(D[1],D[2])), D[7])
ct(mf,F,g) = my(e,pp); e=mfslashexpansion(mf,F,g,0,0,&pp); if(pp[1]!=0, return(0.)); e[1]*1.
nrm(v) = my(s); s=sqrt(sum(i=1,#v, abs(v[i])^2)); if(s==0., v, v/s)
gram3(a,b,c) = my(M,A,B,C); A=nrm(a); B=nrm(b); C=nrm(c); M=matrix(3,3); M[1,1]=A*conj(A~); M[1,2]=A*conj(B~); M[1,3]=A*conj(C~); M[2,1]=conj(M[1,2]); M[2,2]=B*conj(B~); M[2,3]=B*conj(C~); M[3,1]=conj(M[1,3]); M[3,2]=conj(M[2,3]); M[3,3]=C*conj(C~); abs(matdet(M))
gram2(a,b) = my(M,A,B); A=nrm(a); B=nrm(b); M=matrix(2,2); M[1,1]=A*conj(A~); M[1,2]=A*conj(B~); M[2,1]=conj(M[1,2]); M[2,2]=B*conj(B~); abs(matdet(M))
Lv(f,lab,s) = if(f==1, lfun(1,s), lfun(lfuncreate([ZS[f],lab]),s))
pnm(k,u,lab,uodd,uord) = if(u==1, Str("zeta(",k-1,")"), if(uord==2, Str("L(",k-1,",chi_",if(uodd,-u,u),")"), Str("L(",k-1,",chi_",u,"#",lab,")")))
isint(k,D) = if(D[4]!=1, 0, if(k==3, D[3], 1-D[3]))
run(N,k) = my(DD,d,MFS,r1,r2,i,mf,F,g2,pn,keys,Ev,rk,ann,msg); DD=dirs(N,k); d=#DD; if(d==0, return(0)); g2=[0,-1;1,0]; r1=vector(d); r2=vector(d); for(i=1,d, mf=if(N==1, mfinit([N,k],3), mfinit([N,k,Mod(DD[i][8],N)],3)); F=mkE(k,DD[i]); r1[i]=ct(mf,F,[1,0;0,1]); r2[i]=ct(mf,F,g2)); rk=if(gram2(r1,r2)>1e-25, 2, 1); ann=d-rk; print("G1ANN ",N," ",k," ",d," ",ann," rk=",rk); for(i=1,d, print("  ROW ",N," ",k," ",i," psi=[",DD[i][1],",",DD[i][2],",odd",DD[i][3],"] phi=[",DD[i][4],",",DD[i][5],",odd",DD[i][6],"] d=",DD[i][7]," eps=",DD[i][8]," int=",isint(k,DD[i])," a0inf=",if(abs(r1[i])>1e-25,1,0)," a00=",if(abs(r2[i])>1e-25,1,0)," nm=",pnm(k,DD[i][1],DD[i][2],DD[i][3],cord(DD[i][1],DD[i][2])))); keys=Set(vector(d,i,if(isint(k,DD[i]),pnm(k,DD[i][1],DD[i][2],DD[i][3],cord(DD[i][1],DD[i][2])),""))); for(j=1,#keys, if(keys[j]!="", print("  PER ",N," ",k," ",keys[j]," survives=",if(gram3(r1,r2,vector(d,i,if(isint(k,DD[i]) && pnm(k,DD[i][1],DD[i][2],DD[i][3],cord(DD[i][1],DD[i][2]))==keys[j], (-1/2)*DD[i][7]^(1-k)*1., 0.)))>1e-25,1,0)))); Ev=vector(d,i,if(DD[i][4]!=1 && DD[i][6]==1, Lv(DD[i][1],DD[i][2],k-1)*Lv(DD[i][4],DD[i][5],0)*DD[i][7]^(1-k), 0.)); print("  ELEM ",N," ",k," survives=",if(gram3(r1,r2,Ev)>1e-25,1,0)," nonzero=",if(gram2(Ev,Ev)>=0 && sqrt(sum(i=1,d,abs(Ev[i])^2))>1e-25,1,0)); 1
GZ=[1,2,3,4,5,6,7,8,9,10,12,13,16,18,25];
for(t=1,#GZ, N=GZ[t]; for(k=3,4, run(N,k)));
quit;
