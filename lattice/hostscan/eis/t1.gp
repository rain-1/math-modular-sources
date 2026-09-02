f(x) = my(a=1); a+x
print(f(3));
G=znstar(5,1);
print(G.cyc);
print([znconreyconductor(G,n) | n<-[1,2,3,4]]);
print([zncharisodd(G,n) | n<-[1,2,3,4]]);
G1=znstar(1,1);
print(G1.cyc);
print(znconreyconductor(G1,1));
print(zncharisodd(G1,1));
L=lfuncreate([G,2]); print(lfun(L,3));
L2=lfuncreate([G,4]); print(lfun(L2,3));
print(mfdim(mfinit([5,4],3)));
print(mfdim(mfinit([5,3],3)));
quit;
