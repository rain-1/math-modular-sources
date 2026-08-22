
\\ 20_hunt2.gp -- deep identification of the unidentified xi_2 of the aq=0, c=-256 cells.
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_hunt2.gp
default(parisizemax,"14G");
default(realprecision,60);
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");

NN = 700;
PR = 3800;
HT = 100000000000000;

arow(aq,bq,cq,c,r,N) = { my(u=vector(N+2)); u[1]=0; u[2]=1; for(n=0,N-1, u[n+3] = ((aq*n^2+bq*n+cq)*u[n+2] - c*(n-r)^2*u[n+1])/(n+1)^2); vector(N+1,j,u[j+1]); }
brow(aq,bq,cq,c,r,N) = { my(u=vector(N+2)); u[1]=0; for(j=0,min(r,N), u[j+2]=0); if(r+1<=N, u[r+3]=1); for(n=r+1,N-1, u[n+3] = ((aq*n^2+bq*n+cq)*u[n+2] - c*(n-r)^2*u[n+1])/(n+1)^2); vector(N+1,j,u[j+1]); }
xiof(v,N) = { my(aa=arow(v[1],v[2],v[3],v[4],v[5],N), bb=brow(v[1],v[2],v[3],v[4],v[5],N)); bb[N+1]/aa[N+1]; }
small(l) = #l>0 && vecmax(abs(Vec(l)))<HT;

X1 = xiof([0,8,4,-256,0],NN);
X2 = xiof([0,-8,4,-256,2],NN);
X3 = xiof([0,-2568,1284,-256,2],NN);
X4 = xiof([0,8,4,-64,0],NN);
XC = xiof([-32,0,4,256,1],NN);

print("R1 = (0,8,4,-256,0):  (n+1)^2 u_{n+1} = (8n+4) u_n + 256 n^2 u_{n-1}");
print("  A_0..A_11 = ", vector(12,j, arow(0,8,4,-256,0,12)[j]));
print("  B_0..B_7  = ", vector(8,j, brow(0,8,4,-256,0,12)[j]));
print("  xi_2 mod 2^100 = ", lift(X1+O(2^100)));

print("");
print("--- mutual affine relations (lindep, height < 10^14) ---");
{
print("  [1, xi(R1), xi(R2)] : ", if(small(lindep([1+O(2^PR),X1+O(2^PR),X2+O(2^PR)])), lindep([1+O(2^PR),X1+O(2^PR),X2+O(2^PR)]), "none"));
print("  [1, xi(R1), xi(R3)] : ", if(small(lindep([1+O(2^PR),X1+O(2^PR),X3+O(2^PR)])), lindep([1+O(2^PR),X1+O(2^PR),X3+O(2^PR)]), "none"));
print("  [1, xi(R1), xi(R4)] : ", if(small(lindep([1+O(2^PR),X1+O(2^PR),X4+O(2^PR)])), lindep([1+O(2^PR),X1+O(2^PR),X4+O(2^PR)]), "none"));
print("  [1, xi(R1), xi(Cal)]: ", if(small(lindep([1+O(2^PR),X1+O(2^PR),XC+O(2^PR)])), lindep([1+O(2^PR),X1+O(2^PR),XC+O(2^PR)]), "none"));
}

print("");
print("--- algebraicity of xi(R1), degree <= 4 ---");
{
for(d=1,4, my(l=lindep(vector(d+1,j,X1^(j-1)+O(2^PR))));
  print("  degree ",d,": ", if(small(l), l, "none of height < 10^14")));
}

print("");
print("--- three-term with zeta_2(2), zeta_2(3) ---");
{
my(z2=Lp(2,triv,2,PR+20), z3=Lp(2,triv,3,PR+20));
my(l=lindep([1+O(2^PR), z2+O(2^PR), z3+O(2^PR), X1+O(2^PR)]));
print("  [1, zeta_2(2), zeta_2(3), xi(R1)] : ", if(small(l), l, "none of height < 10^14"));
}

print("");
print("--- affine in L_2(s,chi), even quadratic chi of conductor <= 104, s = 2..5 ---");
{
my(hit=0);
for(d=2,104,
  if(!isfundamental(d), next);
  my(ch=[d, vector(d,a, if(gcd(a,d)==1, kronecker(d,a), 0))]);
  for(s=2,5,
    my(lv=Lp(2,ch,s,PR+20));
    if(lv==0, next);
    if(valuation(lv,2)>PR-60, next);
    my(l=lindep([1+O(2^PR), lv+O(2^PR), X1+O(2^PR)]));
    if(small(l) && l[3]!=0,
      my(g=-(l[1]+l[2]*lv)/l[3]);
      if(valuation(X1-g,2)>PR-80,
        print("  HIT: xi(R1) = (",-l[1]," + ",-l[2]," L_2(",s,",kron_",d,"))/",l[3]);
        hit=1))));
my(lv=Lp(2,triv,2,PR+20));
my(l=lindep([1+O(2^PR), lv+O(2^PR), X1+O(2^PR)]));
if(small(l), print("  HIT with zeta_2(2): ", l));
if(!hit, print("  none"));
}

print("");
print("--- 2-adic logarithms of small odd units ---");
{
my(hit=0);
forstep(m=3,31,2,
  my(lg=log(m+O(2^PR)));
  my(l=lindep([1+O(2^PR), lg, X1+O(2^PR)]));
  if(small(l) && l[3]!=0,
    my(g=-(l[1]+l[2]*lg)/l[3]);
    if(valuation(X1-g,2)>PR-80, print("  HIT with log_2(",m,"): ", l); hit=1)));
if(!hit, print("  none"));
}
quit;
