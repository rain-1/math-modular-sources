\\ rows_05_ident.gp -- TARGETED identification of the three unidentified Sym^1 periods
\\   AZ(9,3,-27): 0.1455144820...   Cooper s10: 0.3169253592...   Cooper s18: 0.4842375536...
\\ Small structured bases only:
\\  (a) L(f,s), s=1,2, real+imag parts, for ALL weight-3 newforms with odd nebentypus
\\      at the levels structurally attached to each row;
\\  (b) periods/units: pi^j, sqrt(D), log(fundamental unit), Dirichlet L-values,
\\      zeta(2), zeta(3), Catalan;
\\  (c) two cheap Gamma combinations (Chowla-Selberg for disc -3 and -4).
\\ Tests: x/e rational (2-term), x = r*e + s (3-term), and one joint lindep per basis.
default(parisizemax,8000000000);
PREC=180;
default(realprecision,PREC+40);

\\ ---- recompute the three periods by exact Casoratian telescoping ----
nms=["AZ(9,3,-27)","Cooper s10","Cooper s18"];
P1=[[72,36,6],[24,12,2],[56,28,6]];
P0=[[-432,432,-108],[-256,256,-60],[768,-768,180]];
BB=[6,2,6];
NN=[260,420,2300];
XI=vector(3);
{for(i=1,3,
  my(c=P1[i], d=P0[i], B=BB[i], N=NN[i], a=vector(N+2), W, xi, p1, p0);
  a[1]=1; a[2]=B;
  for(n=1,N, p1=c[1]*n^2+c[2]*n+c[3]; p0=d[1]*n^2+d[2]*n+d[3];
             a[n+2]=(p1*a[n+1]-p0*a[n])/(n+1)^2);
  W=1; xi=0.;
  for(m=0,N, xi += W*1./(a[m+1]*a[m+2]);
             W = W*(d[1]*(m+1)^2+d[2]*(m+1)+d[3])/(m+2)^2);
  XI[i]=xi; print(nms[i]," xi = ",xi));}

\\ ---- weight-3 newform L-values at a given list of levels ----
wt3basis(LS) = { my(C=List());
  for(li=1,#LS, my(L=LS[li], G=znstar(L,1), seen=vectorsmall(L));
    for(m=1,L,
      if(gcd(m,L)!=1 || seen[m], next);
      my(o=znorder(Mod(m,L)));
      for(j=1,o, if(gcd(j,o)==1, seen[lift(Mod(m,L)^j)]=1));
      if(chareval(G, znconreylog(G,m), -1) != 1/2, next);
      iferr(
        my(mf=mfinit([L,3,[G,znconreylog(G,m)]],0));
        if(mfdim(mf)>0, my(EB=mfeigenbasis(mf));
          for(e=1,#EB, my(lf=lfunmf(mf,EB[e]));
            lf = iferr(lfun(lf,2); [lf], E, lf);
            if(type(lf)!="t_VEC", lf=[lf]);
            for(u=1,#lf, for(s=1,2,
              my(v=lfun(lf[u],s), nn=Str("L(f_",L,"_c",m,"_n",e,"_e",u,",",s,")"));
              listput(C,[real(v),nn]);
              if(abs(imag(v))>1e-40, listput(C,[imag(v),Str("Im ",nn)])))))),
      E, )));
  Vec(C); }

constbasis(D1,D2) = { my(C=List());
  listput(C,[Pi,"pi"]); listput(C,[Pi^2,"pi^2"]); listput(C,[Pi^3,"pi^3"]);
  listput(C,[zeta(2),"zeta(2)"]); listput(C,[zeta(3),"zeta(3)"]);
  listput(C,[Catalan,"G"]);
  listput(C,[log(2),"log2"]); listput(C,[log(3),"log3"]);
  listput(C,[log(1+sqrt(2)),"log(1+sqrt2)"]);
  listput(C,[log((1+sqrt(5))/2),"log(golden)"]);
  listput(C,[log(2+sqrt(3)),"log(2+sqrt3)"]);
  listput(C,[sqrt(2),"sqrt2"]); listput(C,[sqrt(3),"sqrt3"]); listput(C,[sqrt(5),"sqrt5"]);
  listput(C,[sqrt(7),"sqrt7"]);
  for(j=1,#[0], );
  my(DS=[-3,-4,-7,-8,-11,5,8,12,-15,-20,-24]);
  for(j=1,#DS, my(D=DS[j], l=lfuncreate(D));
    listput(C,[lfun(l,2),Str("L(2,chi_",D,")")]);
    listput(C,[lfun(l,3),Str("L(3,chi_",D,")")]));
  \\ Chowla-Selberg style Gamma combinations
  listput(C,[gamma(1/3)^6/Pi^4,"Gamma(1/3)^6/pi^4"]);
  listput(C,[gamma(1/4)^4/Pi^3,"Gamma(1/4)^4/pi^3"]);
  listput(C,[gamma(1/8)*gamma(3/8)/(gamma(5/8)*gamma(7/8)),"G(1/8)G(3/8)/G(5/8)G(7/8)"]);
  Vec(C); }

report(x, C, tag) = {
  my(hit=0, tol=10.0^(-PREC+10), bnd=10^14);
  for(j=1,#C, my(e=C[j][1], r);
    if(abs(e)<1e-60, next);
    r=lindep([x,e]);
    if(vecmax(abs(r))<bnd && abs(r[1]*x+r[2]*e)<tol*abs(x),
       print("   HIT (2-term) ",tag,": ",r[1],"*xi + ",r[2],"*",C[j][2]," = 0"); hit=1));
  for(j=1,#C, my(e=C[j][1], r);
    if(abs(e)<1e-60, next);
    r=lindep([x,e,1]);
    if(vecmax(abs(r))<10^10 && abs(r[1]*x+r[2]*e+r[3])<tol*abs(x) && r[2]!=0,
       print("   HIT (3-term) ",tag,": ",r[1],"*xi + ",r[2],"*",C[j][2]," + ",r[3]," = 0"); hit=1));
  hit; }

\\ joint lindep only over a SMALL hand-picked basis (<=6 elements), otherwise meaningless
joint(x, C, tag) = { my(v=concat([x,1.],vector(#C,j,C[j][1])), r);
  if(#C>6, print("   (joint lindep skipped: basis too large to be meaningful)"); return(0));
  r=lindep(v);
  print("   joint lindep over [xi,1,",#C," elts] ",tag,": coeffs = ",r~,
        "  residual = ", abs(sum(j=1,#v, r[j]*v[j])));
  }

LEV=[[3,9,12,27,36,48,72],[5,8,10,16,20,40,80],[3,6,9,12,18,24,36,72]];
\\ hand-picked 4-element joint bases: pi-power, sqrt(disc), log(fund unit), rel. Dirichlet L
SMALL=[[[Pi^2,"pi^2"],[sqrt(3),"sqrt3"],[log(2+sqrt(3)),"log(2+sqrt3)"],[lfun(lfuncreate(-3),2),"L(2,chi-3)"]],[[Pi^2,"pi^2"],[sqrt(5),"sqrt5"],[log((1+sqrt(5))/2),"log(golden)"],[lfun(lfuncreate(-8),2),"L(2,chi-8)"]],[[Pi^2,"pi^2"],[sqrt(2),"sqrt2"],[log(1+sqrt(2)),"log(1+sqrt2)"],[lfun(lfuncreate(-3),2),"L(2,chi-3)"]]];
{for(i=1,3,
  print("");
  print("######## ",nms[i]," ########");
  my(x=XI[i], Cw=wt3basis(LEV[i]), Cc=constbasis(0,0), h1, h2);
  print("  weight-3 L-value basis: ",#Cw," values at levels ",LEV[i]);
  print("  constant basis: ",#Cc," values");
  h1 = report(x,Cw,"wt3");
  h2 = report(x,Cc,"const");
  joint(x, SMALL[i], "small structured");
  if(h1==0 && h2==0, print("   ==> NO relation found at ",PREC," digits (clean negative)")));}
\q
