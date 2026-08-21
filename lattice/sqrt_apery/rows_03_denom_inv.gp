\\ rows_03_denom_inv.gp -- companion denominators (k), invariants, Casoratian
\\ for the eight non-Apery Sym^1 rows.  Exact rational arithmetic throughout.
\\ Row: (n+1)^2 a_{n+1} = p1(n) a_n - p0(n) a_{n-1},  a_0=1, a_1=B ; b_0=0,b_1=1.
default(parisizemax,6000000000);
N=200;
dn(n)=if(n<1,1,lcm(vector(n,i,i)));

nms=["Domb        ","T           ","AZ(9,3,-27) ","AZ(11,5,125)","AZ(7,3,81)  ","Cooper s7   ","Cooper s10  ","Cooper s18  "];
P1=[[20,10,2],[24,12,2],[72,36,6],[88,44,10],[56,28,6],[26,13,2],[24,12,2],[56,28,6]];
P0=[[64,-64,16],[16,-16,4],[-432,432,-108],[2000,-2000,500],[1296,-1296,324],[-27,27,-6],[-256,256,-60],[768,-768,180]];
BB=[2,2,6,10,6,2,2,6];
CHP=[[20,64],[24,16],[72,-432],[88,2000],[56,1296],[26,-27],[24,-256],[56,768]];

{for(i=1,8,
  my(nm=nms[i], c=P1[i], d=P0[i], B=BB[i]);
  my(a=vector(N+2), b=vector(N+2), W=vector(N+1), p1, p0);
  a[1]=1; a[2]=B; b[1]=0; b[2]=1;
  for(n=1,N,
    p1 = c[1]*n^2+c[2]*n+c[3];
    p0 = d[1]*n^2+d[2]*n+d[3];
    a[n+2]=(p1*a[n+1]-p0*a[n])/(n+1)^2;
    b[n+2]=(p1*b[n+1]-p0*b[n])/(n+1)^2);
  W[1]=1;
  for(n=1,N, W[n+1]=W[n]*(d[1]*n^2+d[2]*n+d[3])/(n+1)^2);
  print("=== ",nm," ===");
  print("  a_n n=0..7 : ",vector(8,j,a[j]));
  print("  a_n integral to n=",N,"? ",vector(N+1,j,denominator(a[j]))==vector(N+1,j,1));
  print("  d_n^2 b_n integral n<=",N,"? ",
        vector(N,n,denominator(dn(n)^2*b[n+1]))==vector(N,n,1));
  print("  d_n^1 b_n integral n<=",N,"? ",
        vector(N,n,denominator(dn(n)*b[n+1]))==vector(N,n,1));
  my(f=0); for(n=1,N, if(f==0 && denominator(dn(n)*b[n+1])!=1, f=n));
  print("  first n with d_n*b_n NOT integral: ",f);
  print("  Casoratian == prod_{j<=n} p0(j)/(j+1)^2 for n<=",N,"? ",
        vector(N+1,j,a[j]*b[j+1]-a[j+1]*b[j])==W);
  print("  denominator(b_n) factored, n=1..20:");
  for(n=1,20, print("     n=",n,"  den=",factor(denominator(b[n+1]))));
  print("  d_n^2/den(b_n), n=1..20: ",vector(20,n,dn(n)^2/denominator(b[n+1])));
  my(T=CHP[i][1], D=CHP[i][2], disc=T^2-4*D);
  print("  char poly: x^2 - ",T,"*x + ",D,"    disc = ",disc);
  if(disc>=0,
    my(l1=(T+sqrt(disc))/2, l2=(T-sqrt(disc))/2, tmp);
    if(abs(l2)>abs(l1), tmp=l1; l1=l2; l2=tmp);
    printf("  lambda_1=%.6f  lambda_2=%.6f  c=%d  k=2  score=%.4f  budget=%.4f\n",
           l1,l2,D,log(1/abs(l2))-2, log(abs(l1))-2)
  ,
    printf("  COMPLEX roots, |lambda|=%.6f  c=%d  k=2  score=n/a  budget=%.4f\n",
           sqrt(D),D, log(sqrt(D))-2));
  print("");
);}
\q
