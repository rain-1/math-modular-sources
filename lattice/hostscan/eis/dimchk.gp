default(parisize,600000000);
one(N,k,m) = print("PDIM ",N," ",k," ",m," ",mfdim([N,k,Mod(m,N)],3))
for(N=2,60, for(k=3,4, for(m=1,N, if(gcd(m,N)==1, if(zncharisodd(znstar(N,1),m)==(k%2), one(N,k,m))))));
print("PDIM 1 4 0 ",mfdim([1,4],3));
print("PDIM 1 3 0 ",mfdim([1,3],3));
quit;
