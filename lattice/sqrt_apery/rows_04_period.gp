\\ rows_04_period.gp -- periods xi = lim b_n/a_n of the Sym^1 rows, by EXACT
\\ Casoratian telescoping:  W_n = a_n b_{n+1} - a_{n+1} b_n = prod_{j=1}^n p0(j)/(j+1)^2
\\ (W_0 = 1), hence   xi = sum_{m>=0} W_m/(a_m a_{m+1}).
\\ a_m are exact integers, W_m exact rationals; only the final summation is floating.
\\ Rows with complex characteristic roots (AZ(11,5,125), AZ(7,3,81)) have NO
\\ archimedean limit and are excluded.
default(parisizemax,8000000000);
default(realprecision,320);

nms=["Domb","T","AZ(9,3,-27)","Cooper s7","Cooper s10","Cooper s18"];
P1=[[20,10,2],[24,12,2],[72,36,6],[26,13,2],[24,12,2],[56,28,6]];
P0=[[64,-64,16],[16,-16,4],[-432,432,-108],[-27,27,-6],[-256,256,-60],[768,-768,180]];
BB=[2,2,6,2,2,6];
NN=[420,220,260,220,420,2300];

{for(i=1,6,
  my(nm=nms[i], c=P1[i], d=P0[i], B=BB[i], N=NN[i]);
  my(a=vector(N+2), W, xi, p1, p0);
  a[1]=1; a[2]=B;
  for(n=1,N, p1=c[1]*n^2+c[2]*n+c[3]; p0=d[1]*n^2+d[2]*n+d[3];
             a[n+2]=(p1*a[n+1]-p0*a[n])/(n+1)^2);
  W=1; xi=0.;
  for(m=0,N, xi += W*1./(a[m+1]*a[m+2]);
             W = W*(d[1]*(m+1)^2+d[2]*(m+1)+d[3])/(m+2)^2);
  print(nm,"  N=",N);
  print("  xi = ",xi);
  write("lattice/sqrt_apery/rows_xi.txt", Str(nm,"  ",xi));
);}
print("");
print("=== CM identifications: L(f,2) for the weight-3 CM newforms ===");
{
my(dat=[[12,-3,"eta_2^3 eta_6^3"],[8,-8,"eta_1^2 eta_2 eta_4 eta_8^2"],[7,-7,"eta_1^3 eta_7^3"]]);
for(j=1,3,
  my(L=dat[j][1], D=dat[j][2], nmf=dat[j][3]);
  my(mf=mfinit([L,3,D],0), EB=mfeigenbasis(mf));
  print("level ",L," chi_",D,"  dim S_3^new = ",#EB);
  for(e=1,#EB,
    my(F=EB[e], lf=lfunmf(mf,F), v);
    lf = iferr(lfun(lf,2); [lf], E, lf);
    v=lfun(lf[1],2);
    print("  e=",e,"  q-exp ",mfcoefs(F,10),"   L(f,2) = ",v));
);}
\q
