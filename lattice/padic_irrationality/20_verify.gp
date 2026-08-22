
\\ 20_verify.gp -- exact Phase-2 verification and scoring of the survivors of 20_scan.
\\ Row:  (n+1)^2 u_{n+1} = (aq n^2 + bq n + cq) u_n - c (n-r)^2 u_{n-1}
\\   A-row : A_{-1}=0, A_0=1, recurrence from n=0            (the "denominator" row)
\\   B-row : B_j=0 for j<=r, B_{r+1}=1, recurrence from n=r+1 (the second solution)
\\ Stage 1 : k from n<=60, lambda_1 from the char. equation, provisional S_p.
\\ Stage 2 : for S_p > -0.3, exact integrality to n=300, measured sigma_p and
\\           kappa_p, final S_p / theta_p, and identification of xi_p against
\\           Kubota-Leopoldt L_p(s,chi) (lp.gp of lattice/euler_criterion).
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_verify.gp
default(parisizemax,"12G");
default(realprecision,60);
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");

N1 = 60;
N2 = 300;
NM = 250;
THRESH = -0.3;
SCANOUT = "/home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_scan.out";

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
{ kexp(bb,N) =
  my(kmax=0,d=1);
  for(n=1,N, d=lcm(d,n); my(kk=0);
    while(denominator(bb[n+1]*d^kk)!=1 && kk<25, kk++);
    if(kk>kmax,kmax=kk));
  kmax;
}
{ allint(aa,N) = my(f=1); for(n=0,N, if(denominator(aa[n+1])!=1, f=0; break)); f; }
{ lam1of(aq,c) = my(d=aq^2-4*c); if(d>=0, (abs(aq)+sqrt(d*1.))/2, sqrt(abs(c*1.))); }

CHARS = [triv, chim3, chim4, chi12, chi5,
         [8,[1,0,0,0,0,0,-1,0]],
         [8,[1,0,1,0,-1,0,-1,0]],
         [24,[1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0]]];
CHNAM = ["1","chi_-3","chi_-4","chi_12","chi_5","chi_8","chi_-8","chi_24"];

{ ident(xi, pp, PREC) =
  my(out=List());
  for(s=2,3,
   for(i=1,#CHARS,
    my(lv=Lp(pp, CHARS[i], s, PREC));
    if(lv==0, next);
    my(vl=valuation(lv,pp));
    if(vl>PREC-30, next);
    my(rt=xi/lv, v=valuation(rt,pp));
    if(v<-15, next);
    my(sc=pp^max(0,-v), m=pp^(PREC-25));
    my(g=bestappr(Mod(lift(rt*sc),m)));
    if(g!=0 && (type(g)=="t_FRAC" || type(g)=="t_INT"),
      if(abs(numerator(g))<10^5 && abs(denominator(g))<10^5,
        my(chk=valuation(xi - (g/sc)*lv, pp));
        if(chk>PREC-45,
          listput(out, Str(g,"/",sc," * L_",pp,"(",s,",",CHNAM[i],")   [to ",pp,"^",chk,"]")))))));
  Vec(out);
}

{
my(f=fileopen(SCANOUT), ln, rows=List(), tot=0);
while(ln = filereadstr(f),
  if(#ln==0, next);
  if(Vecsmall(ln)[1]==35, next);
  my(v=[eval(x) | x <- select(y->y!="", strsplit(ln," "))]);
  if(#v>=5, listput(rows, [v[1],v[2],v[3],v[4],v[5]])));
fileclose(f);
rows = Set(rows);
print("distinct survivors from 20_scan: ", #rows);

my(keep=List(), nint=0);
for(i=1,#rows,
  my(v=rows[i], aq=v[1], bq=v[2], cq=v[3], c=v[4], r=v[5]);
  my(aa=arow(aq,bq,cq,c,r,N1));
  if(!allint(aa,N1), next);
  nint++;
  my(bb=brow(aq,bq,cq,c,r,N1), kk=kexp(bb,N1), l1=lam1of(aq,c));
  my(fc=factor(abs(c))[,1]);
  for(j=1,#fc, my(pp=fc[j], sp=valuation(c,pp)*log(pp)-kk-log(l1));
    if(sp>THRESH, listput(keep,[aq,bq,cq,c,r,pp,kk,l1,sp])));
);
print("integral to n=60: ", nint);
print("cells with provisional S_p > ", THRESH, ": ", #keep);
print("");
print("=== stage 2: exact to n=300, measured sigma_p ===");
print("aq  bq  cq  c  r | p | k | lambda_1 | sigma_p | kappa_p | S_p | theta_p");
my(final=List());
for(i=1,#keep,
  my(t=keep[i], aq=t[1],bq=t[2],cq=t[3],c=t[4],r=t[5],pp=t[6]);
  my(aa=arow(aq,bq,cq,c,r,N2));
  if(!allint(aa,N2), printf("%d %d %d %d %d  -- A_n NOT integral to 300, discarded\n",aq,bq,cq,c,r); next);
  my(bb=brow(aq,bq,cq,c,r,N2), kk=kexp(bb,N2), l1=lam1of(aq,c));
  my(xi=bb[N2+1]/aa[N2+1]);
  my(v1=valuation(xi-bb[151]/aa[151],pp), v2=valuation(xi-bb[NM+1]/aa[NM+1],pp));
  my(sg=(v2-v1)*1./(NM-150));
  \\ kappa_p = -lim v_p(A_n)/n, measured as a slope (negative = p-power CONTENT,
  \\ positive = p-power DENOMINATORS).  S_p = (v_p(c)+kappa_p) log p - k - log lambda_1
  \\               = sigma_p log p - k - kappa_p log p - log lambda_1.
  my(kap=-(valuation(aa[NM+1],pp)-valuation(aa[151],pp))*1./(NM-150));
  my(den=kk+kap*log(pp)+log(l1), sp=sg*log(pp)-den, th=sg*log(pp)/den);
  printf("%d %d %d %d %d | p=%d | k=%d | %9.5f | %7.4f | %7.4f | %+8.5f | %8.6f\n",
     aq,bq,cq,c,r,pp,kk,l1,sg,kap,sp,th);
  if(sp>THRESH && sg>0.5, listput(final,[sp,aq,bq,cq,c,r,pp,sg,kk,l1,th]));
);
print("");
print("=== rows with S_p > ", THRESH, " and a genuine p-adic limit ===");
final = vecsort(Vec(final),1,4);
for(i=1,#final, my(t=final[i]);
  printf("S=%+8.5f  theta=%8.6f  P(n)=%d n^2 + %d n + %d,  c=%d, r=%d,  p=%d, sigma=%.3f, k=%d, lambda_1=%.5f\n",
    t[1],t[11],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10]));
print("");
print("=== identification of xi_p ===");
for(i=1,#final, my(t=final[i], aq=t[2],bq=t[3],cq=t[4],c=t[5],r=t[6],pp=t[7],sg=t[8]);
  my(NB=min(1200,ceil(700/sg)+30));
  my(aa=arow(aq,bq,cq,c,r,NB), bb=brow(aq,bq,cq,c,r,NB));
  my(xi=bb[NB+1]/aa[NB+1], PREC=floor(sg*NB)-30);
  printf("(%d,%d,%d; c=%d, r=%d) p=%d : xi_p to %d^%d\n",aq,bq,cq,c,r,pp,pp,PREC);
  my(x=xi+O(pp^PREC), id=ident(x,pp,PREC));
  if(#id==0, print("    no low-height Kubota-Leopoldt match"),
     for(j=1,#id, print("    xi_p = ", id[j])));
  \\ archimedean companion, when |lambda_2| < |lambda_1|
  my(d=aq^2-4*c);
  if(d>0,
    my(l2=(abs(aq)-sqrt(d*1.))/2);
    if(abs(l2)<lam1of(aq,c)-1e-9,
      print("    xi_infty ~ ", bb[N2+1]/aa[N2+1]*1.0)));
);
}
quit;
