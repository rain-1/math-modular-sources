/* 07_orbits.gp -- the cusp-move orbits of Zagier's six rows */
{ crow(aa,bb,dd,N) = my(u=vector(N+2)); u[1]=1;
  for(n=0,N, my(Pn=aa*(n^2+n)+bb, Qn=dd*n^2, prv=if(n==0,0,u[n]));
      u[n+2]=(Pn*u[n+1]-Qn*prv)/(n+1)^2);
  vector(N+1,k,u[k]); }
{ moves(aa,bb,dd) = my(D=aa^2-4*dd, out=[]);
  if(!issquare(D), return("irrational characteristic roots"));
  my(s=sqrtint(D));
  for(e=0,1, my(lam=(aa+(-1)^e*s)/2, mu=aa-lam);
      out=concat(out, [[mu-2*lam, bb-lam, lam^2-lam*mu]]));
  out; }
nms=["A","B","C","D","E","F"]; dat=[[7,2,-8],[9,3,27],[10,3,9],[11,3,-1],[12,4,32],[17,6,72]];
for(j=1,6, print(nms[j], " (a,b,d)=", dat[j], "  cusp-moves -> ", moves(dat[j][1],dat[j][2],dat[j][3])));
print();
print("third placement of E: (a,b,d)=(0,0,-16), u_n =", crow(0,0,-16,14));
print("binom(2m,m)^2 for m=0..7        :", vector(8,m,binomial(2*(m-1),m-1)^2));
print("integral? ", vecmax(vector(15,k,denominator(crow(0,0,-16,14)[k])))==1);
print();
print("E's other placement (-12,-4,32) vs (-1)^n * E: ", vector(8,k,crow(-12,-4,32,10)[k]-(-1)^(k-1)*crow(12,4,32,10)[k]));
quit;
