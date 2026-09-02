\\ shared helpers (mirrors lattice/hostscan/lib.gp, plus cusp-form source tools)
useries(dv, r, nq) = q*prod(t=1, #dv, eta(q^dv[t] + O(q^nq))^r[t]);
Fseries(dv, r, nq) = 1 - sum(t=1, #dv, r[t]*dv[t]*sum(n=1, (nq-1)\dv[t], sigma(n)*q^(dv[t]*n))) + O(q^nq);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
thetaser(cv, r, nq) = sum(m=1, nq-1, (cv[m+1]/m^r)*q^m) + O(q^nq);
denexp(bb, nn) = my(dn=1, km=0); for(n=1, nn, dn=lcm(dn,n); my(bn=bb[n+1]); if(bn!=0, my(de=denominator(bn), kj=0, t=de); while(t>1 && kj<12, kj++; t=t/gcd(t,dn)); if(kj>km, km=kj))); km;
fitrecR(av, ord, dg) = my(nv=(ord+1)*(dg+1), rows=List()); for(n=0, #av-ord-2, my(row=vector(nv)); for(j=0,ord, for(e=0,dg, row[j*(dg+1)+e+1] = n^e*av[n+j+1])); listput(rows,row)); matker(matconcat(Vec(rows)~));
normPR(P) = my(v=[]); for(j=1,#P, v=concat(v,Vec(P[j]))); my(g=content(v)); vector(#P,j,P[j]/g);
primvec(v) = my(d=1); for(i=1,#v, d=lcm(d,denominator(v[i]))); v=v*d; my(g=0); for(i=1,#v, g=gcd(g,v[i])); if(g!=0, v=v/g); v;
minusbasis(N) = my(mf, d, Bb, ai, M, Km, res); mf = mfinit([N,4],1); d = mfdim(mf); if(d==0, return([])); Bb = mfbasis(mf); ai = mfatkininit(mf,N); M = matrix(d,d); for(j=1,d, my(vv=mftobasis(mf, mfatkin(ai, Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Km = matker(M+1); res = List(); for(j=1, matsize(Km)[2], listput(res, mflinear(mf, Km[,j]))); [mf, Vec(res)];
plusbasis(N) = my(mf, d, Bb, ai, M, Kp, res); mf = mfinit([N,4],1); d = mfdim(mf); if(d==0, return([])); Bb = mfbasis(mf); ai = mfatkininit(mf,N); M = matrix(d,d); for(j=1,d, my(vv=mftobasis(mf, mfatkin(ai, Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Kp = matker(M-1); res = List(); for(j=1, matsize(Kp)[2], listput(res, mflinear(mf, Kp[,j]))); [mf, Vec(res)];
