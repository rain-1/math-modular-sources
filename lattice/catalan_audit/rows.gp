\p 300
default(parisize,4000000000);
G = Catalan;
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
MM=40;
ZZ=zud(MM); ZQ=ZZ[1]; ZP=ZZ[2];
print("=== Zudilin row ===");
{
for(m=1,10, my(q=ZQ[m+1],p=ZP[m+1]);
  printf("m=%2d v2(denQ)=%3d oddpartdenQ=%s v2(denP)=%3d oddpartdenP=%s |G-P/Q|=%.4e\n",
    m, valuation(denominator(q),2), denominator(q)>>valuation(denominator(q),2),
    valuation(denominator(p),2), denominator(p)>>valuation(denominator(p),2),
    abs(G-p/q)));
}
e(m) = min(6*m, 4*m+3+logint(2*m-1,2));
ok=1;
{
for(m=1,MM, my(q=ZQ[m+1],p=ZP[m+1],em=e(m),D=lcm(vector(2*m-1,i,i)));
  if(denominator(2^e(m)*ZQ[m+1])!=1, ok=0; print("FAILQ ",m));
  if(denominator(2^e(m)*lcm(vector(2*m-1,i,i))^2*ZP[m+1])!=1, ok=0; print("FAILP ",m)));
}
print("Zudilin integrality m<=",MM,": ",if(ok,"OK","FAIL"));
phi=(1+sqrt(5))/2;
printf("5*log(phi)=%.8f\n",5*log(phi));
{
for(k=1,4, my(m=10*k);
  printf("m=%3d log|Q_m|/m=%.6f  log|Q_m G-P_m|/m=%.6f\n", m, log(abs(ZQ[m+1]))/m, log(abs(ZQ[m+1]*G-ZP[m+1]))/m));
}
write("zud_exact.txt","");
{
for(n=1,13, my(m=3*n, D=lcm(vector(6*n,i,i)), X=2^e(m)*D^2*ZQ[m+1], Y=2^e(m)*D^2*ZP[m+1]);
  write("zud_exact.txt", n, " ", X, " ", Y));
}
\q
