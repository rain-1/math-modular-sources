
\\ 20_hunt.gp -- hard identification of xi_2 for the S_2 > 0 rows that 20_ident.gp
\\ could not match: are they rational?  affinely equivalent to Calegari's constant?
\\ affine in some Kubota-Leopoldt L_2(s,chi)?  Uses p-adic lindep with large height.
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_hunt.gp
default(parisizemax,"12G");
default(realprecision,60);
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");

NN = 700;
{ arow(aq,bq,cq,c,r,N) =
  my(u=vector(N+2)); u[1]=0; u[2]=1;
  for(n=0,N-1, u[n+3] = ((aq*n^2+bq*n+cq)*u[n+2] - c*(n-r)^2*u[n+1])/(n+1)^2);
  vector(N+1,j,u[j+1]);
}
{ brow(aq,bq,cq,c,r,N) =
  my(u=vector(N+2)); u[1]=0;
  for(j=0,min(r,N), u[j+2]=0);
  if(r+1<=N, u[r+3]=1);
  for(n=r+1,N-1, u[n+3] = ((aq*n^2+bq*n+cq)*u[n+2] - c*(n-r)^2*u[n+1])/(n+1)^2);
  vector(N+1,j,u[j+1]);
}
{ xiof(v,N) = my(aa=arow(v[1],v[2],v[3],v[4],v[5],N), bb=brow(v[1],v[2],v[3],v[4],v[5],N));
  bb[N+1]/aa[N+1]; }

\\ even Dirichlet characters of small conductor, as [f,[chi(1..f)]]
{ mkchi(f,d) = [f, vector(f, a, if(gcd(a,f)==1, kronecker(d,a), 0))]; }
CH = List(); CHN = List();
{ listput(CH, triv); listput(CHN, "1"); }
{ for(i=1,20, my(d=[-3,-4,5,8,-8,12,13,-11,17,-7,21,24,-15,28,29,-20,33,37,40,-24][i]);
     listput(CH, mkchi(abs(d)*if(d%4==1,1,1), d)); listput(CHN, Str("kron_",d))); }
CH = Vec(CH); CHN = Vec(CHN);

ROWS = [[0,8,4,-256,0]];
{ ROWS = concat(ROWS, [[0,-8,4,-256,2]]); }
{ ROWS = concat(ROWS, [[0,-2568,1284,-256,2]]); }
{ ROWS = concat(ROWS, [[0,8,4,-64,0]]); }
CAL = [-32,0,4,256,1];

{
my(N=NN, xc=xiof(CAL,N), PR=3900);
print("Calegari control xi_2 (row ",CAL,")");
for(i=1,#ROWS,
  my(v=ROWS[i], x=xiof(v,N));
  print("\n=== P(n)=",v[1],"n^2+",v[2],"n+",v[3],", c=",v[4],", r=",v[5]," ===");
  my(xp=x+O(2^PR));
  \\ 1. rationality, large height
  my(l=lindep([1+O(2^PR), xp]));
  print("  lindep[1,xi]      = ", l, "   (rational iff both entries small)");
  \\ 2. affine relation to Calegari's constant
  my(l2=lindep([1+O(2^PR), xc+O(2^PR), xp]));
  print("  lindep[1,xi_Cal,xi] = ", l2);
  if(#l2==3 && l2[3]!=0 && vecmax(abs(l2))<10^12,
     my(g=-(l2[1]+l2[2]*xc)/l2[3]);
     print("     -> xi = ",-l2[1],"/",l2[3]," + ",-l2[2],"/",l2[3]," * xi_Cal ;  check v_2 = ",
           valuation(x-g,2), " of ", PR));
  \\ 3. affine in a Kubota-Leopoldt value
  my(hit=0);
  for(s=2,3, for(j=1,#CH,
    my(lv=Lp(2,CH[j],s,PR+20));
    if(lv==0, next);
    if(valuation(lv,2)>PR-40, next);
    my(rel=lindep([1+O(2^PR), lv+O(2^PR), xp]));
    if(#rel==3 && rel[3]!=0 && vecmax(abs(rel))<10^12,
      my(g=-(rel[1]+rel[2]*lv)/rel[3]);
      if(valuation(x-g,2)>PR-60,
        print("  xi_2 = (",-rel[1]," + ",-rel[2]," L_2(",s,",",CHN[j],"))/",rel[3],
              "   [v_2 check ",valuation(x-g,2),"]"); hit=1)))); 
  if(!hit, print("  no affine L_2(s,chi) relation, height < 10^12, conductor <= 40, s=2,3"));
);
}
quit;
