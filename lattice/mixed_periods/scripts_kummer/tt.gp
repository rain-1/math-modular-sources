Lt=vector(20);Lt[1]=1;for(n=1,19,Lt[n+1]=lcm(Lt[n],n)); L(n)=if(n<=0,1,Lt[n+1]);
tB(n)=L(n);
f(g)=g(5);
print(tB(5)); print(f(tB)); print(f(n->L(n)));
quit;
