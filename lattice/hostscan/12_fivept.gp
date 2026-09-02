/* The three unit-lambda_2 Fricke configurations that are NOT CDT-shape:
   pole of x at a cusp pair but with extra singular points. Order-4 recurrences. */
default(parisizemax, 12000000000);
default(realprecision, 200);
NQ = 120; NA = 100; NB = 260;
useries(dv, r, nq) = q*prod(t=1, #dv, eta(q^dv[t] + O(q^nq))^r[t]);
Fseries(dv, r, nq) = 1 - sum(t=1, #dv, r[t]*dv[t]*sum(n=1, (nq-1)\dv[t], sigma(n)*q^(dv[t]*n))) + O(q^nq);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
fitrecR(av, ord, dg) = my(nv=(ord+1)*(dg+1), rows=List()); for(n=0, #av-ord-2, my(row=vector(nv)); for(j=0,ord, for(e=0,dg, row[j*(dg+1)+e+1] = n^e*av[n+j+1])); listput(rows,row)); matker(matconcat(Vec(rows)~));
denexp(bb, nn) = my(dn=1, km=0); for(n=1, nn, dn=lcm(dn,n); my(bn=bb[n+1]); if(bn!=0, my(de=denominator(bn), kj=0, t=de); while(t>1 && kj<12, kj++; t=t/gcd(t,dn)); if(kj>km, km=kj))); km;
freetest(a, nn) = my(f=1); for(n=1, nn, if(a[n+1] % (n+1) != 0, f=0; break)); f;
iterrecG(P, ini, nb) = my(o=#P-1, v=vector(nb+o+1)); for(i=1,#ini, v[i]=ini[i]); for(n=0, nb, my(pl=subst(P[o+1],x,n)); if(pl==0, break); my(s=0); for(j=1,o, s += subst(P[j],x,n)*v[n+j]); v[n+o+1] = -s/pl); v;
doit(dv,r,C,B,nm) = my(us=useries(dv,r,NQ), Fs=Fseries(dv,r,NQ), xs=us/(1+B*us+C*us^2)); my(a=peel2(Fs,xs,NA,NQ)); print("== ",nm); print("  a = ",vector(10,i,a[i])); my(o=0,d=0,ker=0); for(oo=2,6, for(dd=3,8, my(kk=fitrecR(a,oo,dd)); if(matsize(kk)[2]==1 && o==0, o=oo; d=dd; ker=kk))); if(o==0, print("  no recurrence found"), print("  minimal recurrence: order ",o," (=",o+1,"-term), coeff degree ",d); my(P=vector(o+1,j,sum(e=0,d,ker[(j-1)*(d+1)+e+1,1]*x^e))); my(g=0); for(j=1,o+1, g=gcd(g,content(P[j]))); P=vector(o+1,j,P[j]/g); my(cp=sum(j=1,o+1, polcoeff(P[j],d)*y^(j-1))); print("  char poly (singularities 1/root): ",cp); print("  roots: ",polroots(cp)); print("  freeint = ",freetest(a,60)); my(ini=vector(o)); ini[1]=0; ini[2]=1; my(bb=iterrecG(P,ini,NB)); my(ini2=vector(o,i,a[i])); my(aa=iterrecG(P,ini2,NB)); my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("  reproduces = ",rep); print("  k = ",denexp(bb,140)); print("  xi  = ",bb[NB]*1.0/aa[NB]); print("  xi' = ",bb[NB-1]*1.0/aa[NB-1]));
doit([1,2,5,10],[-3,1,-1,3],20,9,"N=10 C=20 B=9  lam=9+-4sqrt5 (c=1)");
doit([1,2,3,4,6,12],[-3,2,1,-1,-2,3],12,7,"N=12 C=12 B=7  lam=7+-4sqrt3 (c=1)");
doit([1,2,3,6,9,18],[-2,1,1,-1,-1,2],6,5,"N=18 C=6 B=5  lam=5+-2sqrt6 (c=1)");
print("zeta(2)=",zeta(2)); print("zeta(3)=",zeta(3));
quit;
