/* Mechanism A: the two source integrands lie in one Beukers family; ratio = (xy/(1-xy))^n */
\p 40
G=Catalan;
KZ(x,y,m) = x^(m-1/2)*(1-x)^m*y^m*(1-y)^(m-1/2)/(1-x*y)^(m+1);
KN(x,y,n) = x^(4*n-1/2)*(1-x)^(3*n)*y^(4*n)*(1-y)^(3*n-1/2)/(1-x*y)^(4*n+1);
/* ratio test at random points: KN(n)/KZ(3n) =? (xy/(1-xy))^n */
{
for(n=1,4, my(e=0);
 for(t=1,200, my(x=random(10^12)/10^12*0.999+0.0005, y=random(10^12)/10^12*0.999+0.0005);
   my(r=KN(x,y,n)/KZ(x,y,3*n), s=(x*y/(1-x*y))^n); e=max(e,abs(r/s-1)));
 printf("n=%d  max rel dev of KN(n)/KZ(3n) from (xy/(1-xy))^n : %.3e\n", n, e));
}
/* numeric double integrals, small n, vs the exact rows */
{
II(f) = intnum(x=0,1, intnum(y=0,1, f(x,y)));
for(n=1,2,
  my(m=3*n, JZ=II((x,y)->KZ(x,y,m)), JN=II((x,y)->KN(x,y,n)));
  my(zr=zudrow(n), nr=nestrow(n), D=lcm(vector(6*n,i,i)), E=ee(m));
  my(LZ=zr[1]*G-zr[2], LN=nr[1]*G-nr[2]);
  printf("n=%d m=%d: 4*(-1)^m*JZ*2^e*D^2 = %.12e   X_nG-Y_n = %.12e   ratio=%.12f\n",
     n,m, 4*(-1)^m*JZ*2^E*D^2, LZ, (4*(-1)^m*JZ*2^E*D^2)/LZ);
  printf("        4^{7n}D^2*JN = %.12e   V_nG-U_n = %.12e   ratio=%.12f   sign(LZ)=%d sign(LN)=%d\n",
     4^(7*n)*D^2*JN, LN, 4^(7*n)*D^2*JN/LN, sign(LZ), sign(LN));
);
}
/* range of w=xy/(1-xy) on the open square */
\q
