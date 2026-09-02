Lt=vector(20);Lt[1]=1;for(n=1,19,Lt[n+1]=lcm(Lt[n],n)); L(n)=if(n<=0,1,Lt[n+1]);
tB(n)=L(n); tD(n)=L(n)*L(n\2);
print(type(tB)," ",type(tD)," ",tB(5)," ",tD(5));
quit;
