rad(k)=my(r=1);fordiv(k,d,if(isprime(d),r*=d));r;
print("Hosts with cyclotomic collapse: N a multiple of k*rad(k), and D (=N-1 for minus fam, N+1 for plus fam) a perfect k-th power");
{for(k=3,6, my(K=k*rad(k)); print("  k=",k," (k*rad(k)=",K,"):");
  forstep(N=K,300000,K,
    my(Dm=N-1, Dp=N+1);
    if(ispower(Dm,k,&n1), print("     (1-",N,"x)^{-1/",k,"}  M-family  D=N-1=",Dm,"=",n1,"^",k));
    if(ispower(Dp,k,&n2), print("     (1+",N,"x)^{-1/",k,"}  P-family  D=N+1=",Dp,"=",n2,"^",k));
  ));}
quit;
