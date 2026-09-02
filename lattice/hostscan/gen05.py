import json
hosts=json.load(open('03_units.json'))
hdr = r'''default(parisizemax, 10000000000);
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
'''
body=[]
for h in hosts:
    if h['deg']>2: continue
    N=h['N']; C=h['C']; B=h['B']; r=h['r']
    dv=[d for d in range(1,N+1) if N%d==0]
    body.append(f'''
print("### N=",{N}," C=",{C}," B=",{B}," deg=",{h['deg']});
dv={dv}; r={r}; C={C}; B={B};
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
if(matsize(ker)[2]==1, P = normPR(polsofR(ker[,1], ord, dg)); \\
   print("charpoly = ", sum(j=1,ord+1, polcoeff(P[j],dg)*y^(j-1))); \\
   print("freeint = ", freetest(a, 70)); \\
   if(ord==2, aa=iterrec2(P,a[1],a[2],NB); bb=iterrec2(P,0,1,NB), aa=iterrec3(P,a[1],a[2],a[3],NB); bb=iterrec3(P,0,1,0,NB)); \\
   my(rep=1); for(i=1,min(#a,40), if(aa[i]!=a[i], rep=0; break)); print("reproduces=",rep); \\
   print("k = ", denexp(bb, 150)); \\
   print("xi = ", bb[NB]*1.0/aa[NB]); \\
   print("xiprev = ", bb[NB-1]*1.0/aa[NB-1]));
''')
open('05_rows.gp','w').write(hdr+''.join(body)+'\nquit;\n')
print("hosts written:", sum(1 for h in hosts if h['deg']<=2))
