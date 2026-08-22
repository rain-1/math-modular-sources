\\ mp_primesupport.gp -- slope-prime support of every second/third-order census row.
default(parisizemax,4000000000);
NN=500;
r2(a,b,c,u1,N)={
  my(av=vector(N+2),bv=vector(N+2)); av[1]=1; av[2]=u1; bv[1]=0; bv[2]=1;
  for(n=1,N,
    av[n+2]=((a*n^2+a*n+b)*av[n+1]-c*n^2*av[n])/(n+1)^2;
    bv[n+2]=((a*n^2+a*n+b)*bv[n+1]-c*n^2*bv[n])/(n+1)^2); [av,bv]};
r3(a,b,c,u1,N)={
  my(av=vector(N+2),bv=vector(N+2)); av[1]=1; av[2]=u1; bv[1]=0; bv[2]=1;
  for(n=1,N,
    av[n+2]=((2*n+1)*(a*n^2+a*n+b)*av[n+1]-c*n^3*av[n])/(n+1)^3;
    bv[n+2]=((2*n+1)*(a*n^2+a*n+b)*bv[n+1]-c*n^3*bv[n])/(n+1)^3); [av,bv]};
slp(rr,p,N)={my(t=rr[2][N+1]/rr[1][N+1]-rr[2][N]/rr[1][N]); if(t==0,0,valuation(t,p)*1./(N-1))};
ps=[2,3,5,7,11,13];
show(nm,rr,c)={
  my(sp=vector(#ps,j,slp(rr,ps[j],NN)));
  my(supp=[]); for(j=1,#ps, if(sp[j]>0.5, supp=concat(supp,[ps[j]])));
  print(nm,"  c=",c,"  factor(c)=",if(c!=0,Str(factor(abs(c))),"-"),
        "  sigma_p=",vector(#ps,j,round(sp[j]*1000)/1000.),"  support=",supp);
};
print("=== (R2) Zagier rows: (n+1)^2 u = (a n^2+a n+b)u_n - c n^2 u_{n-1} ===");
show("A (7,2,-8)  ", r2(7,2,-8,2,NN), -8);
show("B (9,3,27)  ", r2(9,3,27,3,NN), 27);
show("C (10,3,9)  ", r2(10,3,9,3,NN), 9);
show("D (11,3,-1) ", r2(11,3,-1,3,NN), -1);
show("E (12,4,32) ", r2(12,4,32,4,NN), 32);
show("F (17,6,72) ", r2(17,6,72,6,NN), 72);
print();
print("=== (R3) third-order rows: (n+1)^3 u = (2n+1)(a n^2+a n+b)u_n - c n^3 u_{n-1} ===");
show("Apery z(3) (17,5,1)  ", r3(17,5,1,5,NN), 1);
show("alpha Domb (10,4,64) ", r3(10,4,64,4,NN), 64);
show("epsilon T (12,4,16)  ", r3(12,4,16,4,NN), 16);
show("zeta AZ(9,3,-27)     ", r3(9,3,-27,3,NN), -27);
show("delta AZ(7,3,81)     ", r3(7,3,81,3,NN), 81);
show("eta AZ(11,5,125)     ", r3(11,5,125,5,NN), 125);
\q
