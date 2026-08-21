\p 300
default(parisize,4000000000);
{
zud(M) = my(Q=vector(M+1), P=vector(M+1), a,b,c);
  Q[1]=1; Q[2]=7/4; P[1]=0; P[2]=13/8;
  for(m=1,M-1,
    a=(2*m+1)^2*(2*m+2)^2*(20*m^2-8*m+1);
    b=3520*m^6+5632*m^5+2064*m^4-384*m^3-156*m^2+16*m+7;
    c=(2*m-1)^2*(2*m)^2*(20*m^2+32*m+13);
    Q[m+2]=(b*Q[m+1]+c*Q[m])/a;
    P[m+2]=(b*P[m+1]+c*P[m])/a);
  [Q,P];
}
MM=40; ZZ=zud(MM); ZQ=ZZ[1]; ZP=ZZ[2];
for(k=1,4, my(m=10*k); print("m=",m," v2(P/Q diff)=",valuation(ZP[m+1]/ZQ[m+1]-ZP[m]/ZQ[m],2)));
print("v2(Q_m) at m=10,20,30,40: ", vector(4,k,valuation(ZQ[10*k+1],2)));
print("log|Q_40|/40 = ",log(abs(ZQ[41]*1.))/40," vs 5log(phi)=",5*log((1+sqrt(5))/2));
\q
