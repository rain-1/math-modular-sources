N=40;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n)); L(n)=Lt[n+1];
intg(s)=my(r=O(x^(N+1)));for(k=0,N-1, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
report(name,s)={my(ok1=1,ok11=1,ok2=1,bad=[]);
 for(n=1,N-1, my(dn=denominator(polcoeff(s,n)));
   if(L(n)%dn!=0, ok1=0);
   if((L(n)*L(n\2))%dn!=0, ok11=0; bad=concat(bad,[n]));
   if((L(n)^2)%dn!=0, ok2=0));
 print(name,": [n]:",ok1," [n][n/2]:",ok11," [n]^2:",ok2, "  fails[n][n/2]: ",if(#bad,vector(min(8,#bad),i,bad[i]),"none"));
}
host(name,HA)={
 print("=== host ",name);
 my(lg=log(1-x+O(x^(N+2))));
 report("HA",HA);
 report("HB=HA*int HA/(1-t)",HA*intg(HA/(1-x)));
 report("HC=HA*int HA log/t",HA*intg(HA*lg/x));
 report("HD=HA*int HA log/(1-t)",HA*intg(HA*lg/(1-x)));
 report("HA*int HA log^2/t",HA*intg(HA*lg^2/x));
 report("HA*int HA log^2/(1-t)",HA*intg(HA*lg^2/(1-x)));
 report("HA*int HA/t*int(HA/(1-t))/HA",HA*intg(HA/x*intg(HA/(1-x))/HA));
}
host("1/sqrt(1-4x)",1/sqrt(1-4*x+O(x^(N+2))));
host("1/sqrt(1-8x)",1/sqrt(1-8*x+O(x^(N+2))));
host("1/sqrt(1-12x)",1/sqrt(1-12*x+O(x^(N+2))));
host("1/sqrt((1-x)(1-9x))",1/sqrt((1-x)*(1-9*x)+O(x^(N+2))));
host("Delannoy 1/sqrt(1-6x+x^2)",1/sqrt(1-6*x+x^2+O(x^(N+2))));
host("trinomial 1/sqrt((1+x)(1-3x))",1/sqrt((1+x)*(1-3*x)+O(x^(N+2))));
