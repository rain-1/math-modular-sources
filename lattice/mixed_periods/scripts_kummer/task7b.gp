default(parisize,"2G");
rad(k)=my(r=1);fordiv(k,d,if(isprime(d),r*=d));r;
ok(k,NN,sg)={my(H=(1-sg*NN*x+O(x^82))^(-1/k));for(j=1,k-1,for(n=0,80,if(denominator(polcoeff(H^j,n))!=1,return(0))));1;}
print("full set of N in 1..400 with (1 -+ N x)^{-j/k} in Z[[x]] for all j=1..k-1  (checked to n=80)");
{for(k=2,7, my(K=k*rad(k));
  forstep(sg=1,-1,-2,
    my(g=List()); for(NN=1,400, if(ok(k,NN,sg), listput(g,NN)));
    my(v=Vec(g), isMult=1); for(i=1,#v, if(v[i]%K!=0, isMult=0));
    my(allMult=1); forstep(NN=K,400,K, if(!setsearch(Set(v),NN), allMult=0));
    print("  k=",k," fam=",if(sg==1,"M","P")," k*rad(k)=",K,": set = ",v, "   == exactly the multiples of k*rad(k)? ", if(isMult&&allMult,"YES","NO"))));}
print();
print("named counterexamples (first non-integral coefficient index of (1 -+ Nx)^{-j/k}, j=1..k-1):");
{my(w=[[3,3,1],[3,3,-1],[3,6,1],[4,4,1],[4,4,-1],[4,2,1],[5,5,1],[6,6,1],[6,12,1],[6,18,1],[7,7,1]]);
 for(i=1,#w, my(k=w[i][1],NN=w[i][2],sg=w[i][3],H=(1-sg*NN*x+O(x^82))^(-1/k));
  print("  (1",if(sg==1,"-","+"),NN,"x)^{-j/",k,"}: ",vector(k-1,j,my(r=-1);for(n=0,80,if(denominator(polcoeff(H^j,n))!=1,r=n;break));r)));}
quit;
