default(parisize,600000000);
cord(f,n) = if(f<=1,1,znorder(Mod(n,f)))
plist(f) = my(G,r); G=znstar(f,1); r=List(); for(n=1,f, if(gcd(n,f)==1, if(type(znconreyconductor(G,n))=="t_INT", listput(r,[n,zncharisodd(G,n),cord(f,n)])))); Vec(r)
PL = vector(60,f,plist(f));
ZS = vector(60,f,znstar(f,1));
indlab(u,lab,N) = if(u==N, lab, znconreyexp(ZS[N], zncharinduce(ZS[u],lab,N)))
pclass(k,uodd) = if(k==3, if(uodd,"I","E"), if(uodd,"E","I"))
pnm(k,u,lab,uodd,uord) = if(u==1, Str("zeta(",k-1,")"), if(uord==2, Str("L(",k-1,",chi_",if(uodd,-u,u),")"), Str("L(",k-1,",chi_",u,"#",lab,")")))
emit(N,k,u,ps,v,ph,d) = my(pl,ql,el,cls,nm); pl=indlab(u,ps[1],N); ql=indlab(v,ph[1],N); el=lift(Mod(pl*ql,N)); cls=if(ph[2]||v==1, pclass(k,ps[2]), "Z"); nm=if(cls=="Z","0",pnm(k,u,ps[1],ps[2],ps[3])); print("DIR ",N," ",k," ",u," ",ps[1]," ",ps[2]," ",ps[3]," ",v," ",ph[1]," ",ph[2]," ",ph[3]," ",d," ",el," ",cls," ",nm)
for(N=1,60, for(k=3,4, par=if(k==3,1,0); fordiv(N,u, for(i=1,#PL[u], ps=PL[u][i]; fordiv(N/u,v, for(j=1,#PL[v], ph=PL[v][j]; if((ps[2]+ph[2])%2==par, fordiv(N/(u*v),d, emit(N,k,u,ps,v,ph,d)))))))));
quit;
