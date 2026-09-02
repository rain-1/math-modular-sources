Lt=vector(20);Lt[1]=1;for(n=1,19,Lt[n+1]=lcm(Lt[n],n)); L(n)=if(n<=0,1,Lt[n+1]);
tB(n)=L(n);
exc(s,typ)={my(mx=1);for(n=0,3,my(dn=denominator(polcoeff(s,n)),T=typ(n),e=dn/gcd(dn,T));1);mx;}
exc2(s,tf)={my(mx=1);for(n=0,3,my(dn=denominator(polcoeff(s,n)));my(T=tf(n));my(e=dn/gcd(dn,T)));mx;}
s=1/2+x/3+O(x^5);
print(exc2(s,tB));
print(exc(s,tB));
quit;
