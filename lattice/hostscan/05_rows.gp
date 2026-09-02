default(parisizemax, 10000000000);
default(realprecision, 210);
NQ = 150;
NA = 116;
NB = 320;
useries(dv, r, nq) = q*prod(t=1, #dv, eta(q^dv[t] + O(q^nq))^r[t]);
Fseries(dv, r, nq) = 1 - sum(t=1, #dv, r[t]*dv[t]*sum(n=1, (nq-1)\dv[t], sigma(n)*q^(dv[t]*n))) + O(q^nq);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
fitrecR(av, ord, dg) = my(nv=(ord+1)*(dg+1), rows=List()); for(n=0, #av-ord-2, my(row=vector(nv)); for(j=0,ord, for(e=0,dg, row[j*(dg+1)+e+1] = n^e*av[n+j+1])); listput(rows,row)); matker(matconcat(Vec(rows)~));
polsofR(kv, ord, dg) = vector(ord+1, j, sum(e=0, dg, kv[(j-1)*(dg+1)+e+1]*x^e));
normPR(P) = my(v=[]); for(j=1,#P, v=concat(v,Vec(P[j]))); my(g=content(v)); vector(#P,j,P[j]/g);
freetest(a, nn) = my(f=1); for(n=1, nn, if(a[n+1] % (n+1) != 0, f=0; break)); f;
denexp(bb, nn) = my(dn=1, km=0); for(n=1, nn, dn=lcm(dn,n); my(bn=bb[n+1]); if(bn!=0, my(de=denominator(bn), kj=0, t=de); while(t>1 && kj<12, kj++; t=t/gcd(t,dn)); if(kj>km, km=kj))); km;
iterrec2(P, a0, a1, nb) = my(v=vector(nb+3)); v[1]=a0; v[2]=a1; for(n=0, nb, my(p0=subst(P[1],x,n), p1=subst(P[2],x,n), p2=subst(P[3],x,n)); if(p2==0, break); v[n+3] = -(p0*v[n+1]+p1*v[n+2])/p2); v;
iterrec3(P, a0, a1, a2, nb) = my(v=vector(nb+4)); v[1]=a0; v[2]=a1; v[3]=a2; for(n=0, nb, my(p0=subst(P[1],x,n), p1=subst(P[2],x,n), p2=subst(P[3],x,n), p3=subst(P[4],x,n)); if(p3==0, break); v[n+4] = -(p0*v[n+1]+p1*v[n+2]+p2*v[n+3])/p3); v;

print("### N=",2," C=",4096," B=",127," deg=",1);
dv=[1, 2]; r=[-24, 24]; C=4096; B=127;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",2," C=",4096," B=",129," deg=",1);
dv=[1, 2]; r=[-24, 24]; C=4096; B=129;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",3," C=",729," B=",53," deg=",1);
dv=[1, 3]; r=[-12, 12]; C=729; B=53;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",3," C=",729," B=",55," deg=",1);
dv=[1, 3]; r=[-12, 12]; C=729; B=55;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",4," C=",256," B=",31," deg=",1);
dv=[1, 2, 4]; r=[-8, 0, 8]; C=256; B=31;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",4," C=",256," B=",33," deg=",1);
dv=[1, 2, 4]; r=[-8, 0, 8]; C=256; B=33;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",6," C=",64," B=",15," deg=",2);
dv=[1, 2, 3, 6]; r=[-6, 6, -6, 6]; C=64; B=15;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",6," C=",64," B=",17," deg=",2);
dv=[1, 2, 3, 6]; r=[-6, 6, -6, 6]; C=64; B=17;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",6," C=",72," B=",17," deg=",1);
dv=[1, 2, 3, 6]; r=[-5, 1, -1, 5]; C=72; B=17;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",6," C=",81," B=",17," deg=",2);
dv=[1, 2, 3, 6]; r=[-4, -4, 4, 4]; C=81; B=17;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",6," C=",81," B=",19," deg=",2);
dv=[1, 2, 3, 6]; r=[-4, -4, 4, 4]; C=81; B=19;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",7," C=",49," B=",13," deg=",1);
dv=[1, 7]; r=[-4, 4]; C=49; B=13;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",7," C=",49," B=",15," deg=",1);
dv=[1, 7]; r=[-4, 4]; C=49; B=15;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",8," C=",16," B=",7," deg=",2);
dv=[1, 2, 4, 8]; r=[-8, 16, -16, 8]; C=16; B=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",8," C=",16," B=",9," deg=",2);
dv=[1, 2, 4, 8]; r=[-8, 16, -16, 8]; C=16; B=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",8," C=",64," B=",15," deg=",2);
dv=[1, 2, 4, 8]; r=[0, -12, 12, 0]; C=64; B=15;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",8," C=",64," B=",17," deg=",2);
dv=[1, 2, 4, 8]; r=[0, -12, 12, 0]; C=64; B=17;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",10," C=",16," B=",7," deg=",2);
dv=[1, 2, 5, 10]; r=[-4, 4, -4, 4]; C=16; B=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",10," C=",16," B=",9," deg=",2);
dv=[1, 2, 5, 10]; r=[-4, 4, -4, 4]; C=16; B=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",10," C=",20," B=",9," deg=",1);
dv=[1, 2, 5, 10]; r=[-3, 1, -1, 3]; C=20; B=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",10," C=",25," B=",9," deg=",2);
dv=[1, 2, 5, 10]; r=[-2, -2, 2, 2]; C=25; B=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",10," C=",25," B=",11," deg=",2);
dv=[1, 2, 5, 10]; r=[-2, -2, 2, 2]; C=25; B=11;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",4," B=",3," deg=",2);
dv=[1, 2, 3, 4, 6, 12]; r=[-7, 12, 5, -5, -12, 7]; C=4; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",4," B=",5," deg=",2);
dv=[1, 2, 3, 4, 6, 12]; r=[-7, 12, 5, -5, -12, 7]; C=4; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",9," B=",5," deg=",2);
dv=[1, 2, 3, 4, 6, 12]; r=[-4, 4, 4, -4, -4, 4]; C=9; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",9," B=",7," deg=",2);
dv=[1, 2, 3, 4, 6, 12]; r=[-4, 4, 4, -4, -4, 4]; C=9; B=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",12," B=",7," deg=",1);
dv=[1, 2, 3, 4, 6, 12]; r=[-3, 2, 1, -1, -2, 3]; C=12; B=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",16," B=",7," deg=",2);
dv=[1, 2, 3, 4, 6, 12]; r=[-2, 0, -2, 2, 0, 2]; C=16; B=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",16," B=",9," deg=",2);
dv=[1, 2, 3, 4, 6, 12]; r=[-2, 0, -2, 2, 0, 2]; C=16; B=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",36," B=",11," deg=",2);
dv=[1, 2, 3, 4, 6, 12]; r=[1, -8, -3, 3, 8, -1]; C=36; B=11;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",12," C=",36," B=",13," deg=",2);
dv=[1, 2, 3, 4, 6, 12]; r=[1, -8, -3, 3, 8, -1]; C=36; B=13;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",15," C=",9," B=",5," deg=",2);
dv=[1, 3, 5, 15]; r=[-2, 2, -2, 2]; C=9; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",15," C=",9," B=",7," deg=",2);
dv=[1, 3, 5, 15]; r=[-2, 2, -2, 2]; C=9; B=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",16," C=",4," B=",3," deg=",2);
dv=[1, 2, 4, 8, 16]; r=[-4, 6, 0, -6, 4]; C=4; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",16," C=",4," B=",5," deg=",2);
dv=[1, 2, 4, 8, 16]; r=[-4, 6, 0, -6, 4]; C=4; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",16," C=",16," B=",7," deg=",2);
dv=[1, 2, 4, 8, 16]; r=[0, -4, 0, 4, 0]; C=16; B=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",16," C=",16," B=",9," deg=",2);
dv=[1, 2, 4, 8, 16]; r=[0, -4, 0, 4, 0]; C=16; B=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",18," C=",4," B=",3," deg=",2);
dv=[1, 2, 3, 6, 9, 18]; r=[-3, 3, 2, -2, -3, 3]; C=4; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",18," C=",4," B=",5," deg=",2);
dv=[1, 2, 3, 6, 9, 18]; r=[-3, 3, 2, -2, -3, 3]; C=4; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",18," C=",6," B=",5," deg=",1);
dv=[1, 2, 3, 6, 9, 18]; r=[-2, 1, 1, -1, -1, 2]; C=6; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",18," C=",9," B=",5," deg=",2);
dv=[1, 2, 3, 6, 9, 18]; r=[-1, -1, 0, 0, 1, 1]; C=9; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",18," C=",9," B=",7," deg=",2);
dv=[1, 2, 3, 6, 9, 18]; r=[-1, -1, 0, 0, 1, 1]; C=9; B=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",22," C=",4," B=",3," deg=",2);
dv=[1, 2, 11, 22]; r=[-2, 2, -2, 2]; C=4; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",22," C=",4," B=",5," deg=",2);
dv=[1, 2, 11, 22]; r=[-2, 2, -2, 2]; C=4; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",24," C=",2," B=",3," deg=",2);
dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[-3, 4, 1, 0, 0, -1, -4, 3]; C=2; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",24," C=",6," B=",5," deg=",2);
dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[-1, 1, -1, -3, 3, 1, -1, 1]; C=6; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",28," C=",4," B=",3," deg=",2);
dv=[1, 2, 4, 7, 14, 28]; r=[-1, 0, 1, -1, 0, 1]; C=4; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",28," C=",4," B=",5," deg=",2);
dv=[1, 2, 4, 7, 14, 28]; r=[-1, 0, 1, -1, 0, 1]; C=4; B=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",30," C=",2," B=",3," deg=",2);
dv=[1, 2, 3, 5, 6, 10, 15, 30]; r=[-2, 2, 1, 1, -1, -1, -2, 2]; C=2; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",32," C=",2," B=",3," deg=",2);
dv=[1, 2, 4, 8, 16, 32]; r=[-2, 3, -1, 1, -3, 2]; C=2; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",40," C=",2," B=",3," deg=",2);
dv=[1, 2, 4, 5, 8, 10, 20, 40]; r=[-1, 0, 2, 1, -1, -2, 0, 1]; C=2; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",46," C=",2," B=",3," deg=",2);
dv=[1, 2, 23, 46]; r=[-1, 1, -1, 1]; C=2; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

print("### N=",48," C=",2," B=",3," deg=",2);
dv=[1, 2, 3, 4, 6, 8, 12, 16, 24, 48]; r=[-1, 1, 1, -1, -2, 2, 1, -1, -1, 1]; C=2; B=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ); xs = us/(1+B*us+C*us^2);
a = peel2(Fs, xs, NA, NQ);
print("a = ", vector(9,i,a[i]));
ok=1; for(t=1,#a, if(denominator(a[t])!=1, ok=0; break)); print("integral=",ok);
ker = fitrecR(a, 2, 3); ord=2; dg=3;
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,4); dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,2,5); dg=5);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,3); ord=3; dg=3);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,4); ord=3; dg=4);
if(matsize(ker)[2]!=1, ker=fitrecR(a,3,5); ord=3; dg=5);
print("record=",ord," recdeg=",dg," recdim=",matsize(ker)[2]);
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \
   print("freeint = ", freetest(a, 70)); \
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \
   print("k = ", denexp(bb, 150)); \
   print("xi = ", bb[NB]*1.0/aa[NB]); \
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));

quit;
