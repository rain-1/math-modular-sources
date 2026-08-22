
\\ 20_analyse.gp -- full diagnostics for the rows that 20_verify.gp scores S_p > 0.
\\ For each row it decides whether the cell is GENUINE or DEGENERATE, by testing
\\   (a) the TRUE archimedean rate lim log max(|A_n|,|B_n|)/n  (not just the char. root),
\\   (b) kappa_p and max v_p(A_n),
\\   (c) nonvanishing: is A_n xi_p - B_n = 0 for large n  (Calegari Lemma 3.2),
\\   (d) is xi_p rational (bestappr in Z_p),
\\   (e) identification of xi_p against Kubota-Leopoldt L_p(s,chi).
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_analyse.gp
default(parisizemax,"12G");
default(realprecision,60);
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");

NN = 600;

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

CHV = [triv, chim3, chim4, chi12, chi5];
{ CHV = concat(CHV, [[8,[1,0,0,0,0,0,-1,0]]]); }
{ CHV = concat(CHV, [[8,[1,0,1,0,-1,0,-1,0]]]); }
{ CHV = concat(CHV, [[24,[1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0]]]); }
CHN = ["1","chi_-3","chi_-4","chi_12","chi_5","chi_8","chi_-8","chi_24"];

{ ident(xi, pp, PREC) =
  my(out=List());
  for(s=2,3,
   for(i=1,#CHV,
    my(lv=Lp(pp, CHV[i], s, PREC));
    if(lv==0, next);
    if(valuation(lv,pp)>PREC-30, next);
    my(rt=xi/lv, v=valuation(rt,pp));
    if(v<-15, next);
    my(sc=pp^max(0,-v), m=pp^(PREC-25));
    my(g=bestappr(Mod(lift(rt*sc),m)));
    if(g!=0,
      if(abs(numerator(g))<10^5 && abs(denominator(g))<10^5,
        my(chk=valuation(xi-(g/sc)*lv,pp));
        if(chk>PREC-45,
          listput(out, Str(g,"/",sc," * L_",pp,"(",s,",",CHN[i],")   [to ",pp,"^",chk,"]")))))));
  Vec(out);
}

{ analyse(aq,bq,cq,c,r,tag) =
  my(N=NN, aa=arow(aq,bq,cq,c,r,N), bb=brow(aq,bq,cq,c,r,N));
  print("\n---- ", tag, " :  P(n) = ",aq,"n^2 + ",bq,"n + ",cq,",  c = ",c,",  r = ",r," ----");
  if(!allint(aa,N), print("  A_n NOT integral to ",N," -- REJECTED"); return(0));
  my(kk=kexp(bb,N));
  my(rt=vector(3,j, my(n=[400,500,600][j]);
       log(max(abs(aa[n+1]*1.),abs(bb[n+1]*1.)))/n));
  print("  A_n integral to ",N,";  k = ",kk,
        ";  measured log max(|A|,|B|)/n at n=400,500,600 = ", rt);
  my(d=aq^2-4*c, l1=if(d>=0,(abs(aq)+sqrt(d*1.))/2,sqrt(abs(c*1.))));
  print("  char. roots of x^2-",aq,"x+",c,":  lambda_1 = ",l1,"  (log = ",log(l1),")");
  my(fc=factor(abs(c))[,1]);
  for(j=1,#fc, my(pp=fc[j]);
    my(xi=bb[N+1]/aa[N+1]);
    my(v1=valuation(xi-bb[301]/aa[301],pp), v2=valuation(xi-bb[501]/aa[501],pp));
    my(sg=(v2-v1)*1./200, kap=-valuation(aa[N+1],pp)*1./N);
    print("  p = ",pp,":  sigma_p = ",sg,"   kappa_p = ",kap,
          "   max v_p(A_n), n<=600 = ", vecmax(vector(N+1,i,if(aa[i]!=0,valuation(aa[i],pp),0))));
    \\ (c) nonvanishing: is B_n/A_n eventually constant (over Q)?
    my(e1=bb[401]/aa[401], e2=bb[501]/aa[501], e3=bb[601]/aa[601]);
    if(e1==e2 && e2==e3,
      print("  *** DEGENERATE: B_n/A_n is exactly constant for n>=400 (A_n xi - B_n = 0);"),
      print("  nonvanishing OK: B_n/A_n differs at n=400,500,600"));
    my(PREC=floor(sg*N)-40);
    if(PREC>60,
      my(x=xi+O(pp^PREC));
      my(g=bestappr(Mod(lift(x*pp^max(0,-valuation(x,pp))), pp^(PREC-25))));
      my(sc=pp^max(0,-valuation(x,pp)));
      if(g!=0 && abs(numerator(g))<10^6 && abs(denominator(g))<10^6 &&
         valuation(xi-g/sc,pp)>PREC-45,
         print("  *** xi_",pp," = ",g/sc," is RATIONAL -- degenerate cell"),
         my(id=ident(x,pp,PREC));
         if(#id==0, print("  xi_p: not rational (height<10^6) and no low-height KL match"),
            for(t=1,#id, print("  xi_",pp," = ",id[t])))));
    my(den=kk+max(kap,0)*log(pp)+log(l1));
    printf("  ==> S_p = %+9.5f,  theta_p = %9.6f\n", sg*log(pp)-den, sg*log(pp)/den);
  );
  1;
}

analyse(-32,0,4,256,1,"CALEGARI's X_1(4) Catalan row (control)");
analyse(-32,32,4,256,2,"aq=-32 family, r=2");
analyse(-32,32,20,256,2,"aq=-32 family, r=2, cq=20");
analyse(-32,64,20,256,3,"aq=-32 family, r=3");
analyse(32,-32,4,256,1,"aq=+32 family, r=1");
analyse(32,0,4,256,0,"aq=+32 family, r=0");
analyse(0,-8,4,-64,2,"aq=0, c=-64 family, small member");
analyse(0,-16,8,-64,2,"aq=0, c=-64 family");
analyse(0,-2568,1284,-64,2,"aq=0, c=-64 family, large member");
analyse(0,-8,4,-256,2,"aq=0, c=-256 family, small member");
analyse(0,-2568,1284,-256,2,"aq=0, c=-256 family, large member");
analyse(0,8,4,-64,0,"aq=0, c=-64, r=0");
analyse(0,8,4,-256,0,"aq=0, c=-256, r=0");
quit;
