
\\ 20_ident.gp -- identify xi_2 for every row that 20_verify.gp scores S_2 > 0.
\\ xi_p is only defined up to the choice of second solution, so it is an AFFINE
\\ function of the Kubota-Leopoldt value, not a multiple: we run p-adic lindep on
\\ [1, L_p(s,chi), xi_p].  A relation with z != 0 gives xi_p = -(x + y*L)/z.
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_ident.gp
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

CHV = [triv, chim3, chim4, chi12, chi5];
{ CHV = concat(CHV, [[8,[1,0,0,0,0,0,-1,0]]]); }
{ CHV = concat(CHV, [[8,[1,0,1,0,-1,0,-1,0]]]); }
{ CHV = concat(CHV, [[24,[1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0]]]); }
CHN = ["1","chi_-3","chi_-4","chi_12","chi_5","chi_8","chi_-8","chi_24"];

{ ROWS = [[-32,0,4,256,1],[-32,32,4,256,2],[-32,32,20,256,2],[-32,64,20,256,3],
        [-32,64,36,256,3],[0,8,4,-256,0],[0,16,8,-256,0],[0,-8,4,-256,2],
        [0,-2568,1284,-256,2],[0,8,4,-64,0],[0,16,8,-64,0],[0,-8,4,-64,2],
        [0,-2568,1284,-64,2],[32,32,12,256,0],[32,0,12,256,1],[32,-32,44,256,2],
        [32,-64,28,256,3]]; }

{
for(i=1,#ROWS,
  my(v=ROWS[i], aq=v[1],bq=v[2],cq=v[3],c=v[4],r=v[5], N=NN);
  my(aa=arow(aq,bq,cq,c,r,N), bb=brow(aq,bq,cq,c,r,N));
  print("\n=== P(n)=",aq,"n^2+",bq,"n+",cq,",  c=",c,",  r=",r," ===");
  if(aa[N+1]==0, print("  A_N = 0, skipped"); next);
  my(xi=bb[N+1]/aa[N+1], pp=2);
  my(v1=valuation(xi-bb[301]/aa[301],pp), v2=valuation(xi-bb[501]/aa[501],pp));
  my(sg=(v2-v1)/200.);
  my(PR=floor(sg*500)-40);
  print("  sigma_2 = ",sg,",  xi_2 known to 2^",PR);
  my(x=xi+O(pp^PR));
  \\ rationality
  my(lr=lindep([1+O(pp^PR), x]));
  if(#lr==2 && lr[2]!=0 && abs(lr[1])<10^8 && abs(lr[2])<10^8,
     my(g=-lr[1]/lr[2]);
     if(valuation(xi-g,pp)>PR-30, print("  *** xi_2 = ",g," RATIONAL -> degenerate cell"); next));
  my(found=0);
  for(s=2,3, for(j=1,#CHV,
    my(lv=Lp(pp,CHV[j],s,PR+10));
    if(lv==0, next);
    my(L=lv+O(pp^PR));
    my(rel=lindep([1+O(pp^PR), L, x]));
    if(#rel==3 && rel[3]!=0,
      my(h=vecmax(abs(rel)));
      if(h<10^6,
        my(g=-(rel[1]+rel[2]*lv)/rel[3]);
        if(valuation(xi-g,pp)>PR-40,
          print("  xi_2 = (",-rel[1]," + ",-rel[2]," * L_2(",s,",",CHN[j],")) / ",rel[3]);
          found=1; break(2))))));
  if(!found, print("  no affine Kubota-Leopoldt relation of height < 10^6 found"));
);
}
quit;
