default(realprecision,130);
M=520;
et(k)=prod(n=1,M+1,1-q^(k*n)+O(q^(M+2)));
h=et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/q;
X=h/((h+3)*(h+4));
D(s)=q*deriv(s,q);
qc=exp(-Pi/sqrt(3));
ev(s)=subst(truncate(s),q,qc);
G13=gamma(1/3); W2=G13^6/Pi^4;
Dh=ev(D(h)); a=Dh/W2;
print("Dh(t*)/W^2 = ",a);
for(d=2,16, my(p=algdep(a,d)); if(vecmax(abs(Vec(p)))<10^12, print("deg ",d,": ",p)));
for(k=[2,3,4,6,12], print("a^",k," algdep deg 4: ",algdep(a^k,4)));
print("a^12 algdep deg 2: ",algdep(a^12,2)," deg 6: ",algdep(a^12,6));
\\ also try a*sqrt(3), a*2^(1/3), a*3^(1/4)
print("(a*3^(1/4))^4 deg 2: ",algdep((a*3^(1/4))^4,2),"  (a*2^(1/3))^3 deg2: ",algdep((a*2^(1/3))^3,2),"  (a*sqrt3)^2 deg 2: ",algdep((a*sqrt(3))^2,2));
print("(a^2*3^(1/2))^2 deg2: ",algdep((a^2*sqrt(3))^2,2));
quit;
