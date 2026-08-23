/* 09_foldreg.gp -- fold-regular subspace and period identification.
   r_n(Phi) = b_n(Phi)/a_n is LINEAR in Phi.  For a fold-regular class the
   convergence r_n -> xi is geometric (rate lambda_2/lambda_1); for a class with
   a log^2 at the fold (or a non-zero constant term) it is O(1/n).  So the
   "slow" directions are spanned by the difference vectors w(n2)-w(n1), and the
   fold-regular subspace is their kernel.  On it we read off xi and run lindep
   against {1, zeta(2), G}. */
default(parisizemax, 8000000000);
default(realprecision, 80);
read("lattice/catalan_al_hosts/lib2.gp");

NQ = 260;
aldiv(n,qq,d) = (qq/gcd(d,qq))*gcd(d,n/qq);
{allinv(qs) = my(r=#qs, res=List(1)); for(mask=1, 2^r-1, my(q=1); for(j=1,r, if(bittest(mask,j-1), q*=qs[j])); listput(res,q)); Vec(res);}
ncu(n,c) = eulerphi(gcd(c,n/c));
{ratproj(mm, ev) = my(d=matsize(mm)[1]);
  my(others = select(z->abs(z-ev)>1e-10, [1,-1,I,-I]));
  my(p1 = matid(d)); for(k=1,#others, p1 = p1*(mm - others[k]*matid(d))/(ev-others[k]));
  matrix(d,d,i,j, bestappr(real(p1[i,j]),10^6));}
{qexpv(mfx, cv, nn) = my(bb=mfbasis(mfx), s=vector(nn+1));
  for(i=1,#bb, if(cv[i]!=0, s += cv[i]*mfcoefs(bb[i],nn))); s;}
vec2ser(v) = sum(n=1,#v, v[n]*x^(n-1)) + O(x^#v);
coefv(s,nn) = vector(nn+1, i, polcoeff(s,i-1));
Dop(s) = my(n=serprec(s,x)); sum(k=1,n-1, k*polcoeff(s,k)*x^k) + O(x^n);
Dm2(s) = my(n=serprec(s,x)); sum(k=1,n-1, polcoeff(s,k)/k^2*x^k) + O(x^n);
dnn(n) = lcm(vector(max(n,1),j,j));
GG = Catalan; Z2 = Pi^2/6;

{analyse(nn, qs, evs, tag) =
  my(dv = divisors(nn), kk = #dv, invs = allinv(qs));
  my(orb = vector(kk,i,i));
  for(i=1,kk, for(j=1,#invs, my(d2 = aldiv(nn,invs[j],dv[i]));
    for(l=1,kk, if(dv[l]==d2, my(m=min(orb[i],orb[l]), o1=orb[i], o2=orb[l]);
      for(z=1,kk, if(orb[z]==o1||orb[z]==o2, orb[z]=m))))));
  my(oinf = orb[kk], labs = Set(Vec(orb)));
  my(mf1 = mfinit([nn,1,-4],4), mf3e = mfinit([nn,3,-4],3), mf3 = mfinit([nn,3,-4],4));
  my(d1 = mfdim(mf1), d3 = mfdim(mf3), d3e = mfdim(mf3e));
  my(p1 = matid(d1), p3 = matid(d3), p3e = matid(d3e));
  for(j=1,#qs,
    my(al1=mfatkininit(mf1,qs[j]), al3=mfatkininit(mf3,qs[j]), al3e=mfatkininit(mf3e,qs[j]));
    p1 = p1*ratproj(al1[2]/al1[3], evs[j]);
    p3 = p3*ratproj(al3[2]/al3[3], evs[j]);
    p3e = p3e*ratproj(al3e[2]/al3e[3], evs[j]));
  my(V1 = matimage(p1), V3e = matimage(p3e));
  my(b3 = mfbasis(mf3), mat3 = matrix(NQ+1, d3, i,j, mfcoefs(b3[j],NQ)[i]));
  for(li=1,#labs, my(lb=labs[li]); if(lb != oinf,
    my(ninf = sum(i=1,kk, if(orb[i]==oinf, ncu(nn,dv[i]), 0)));
    my(npol = sum(i=1,kk, if(orb[i]==lb,   ncu(nn,dv[i]), 0)));
    my(einf = #invs/ninf, epol = #invs/npol);
    if(denominator(einf)!=1 || denominator(epol)!=1, next);
    my(ordv = vector(kk));
    for(i=1,kk, if(orb[i]==oinf, ordv[i]=einf); if(orb[i]==lb, ordv[i]=-epol));
    my(rr = unitexp(nn, ordv));
    if(sum(i=1,kk,rr[i]) != 0, next);
    my(tt = hauptmod(nn, rr, NQ+1), dt); dt = Dop(tt);
    my(nv = matsize(V1)[2], cand = List());
    for(i=1,nv, listput(cand, coefv(vec2ser(qexpv(mf1,V1[,i]~,NQ))*dt + O(x^(NQ+1)), NQ)));
    my(bigm = matrix(NQ+1, nv+d3));
    for(i=1,nv, for(r=1,NQ+1, bigm[r,i] = cand[i][r]));
    for(j=1,d3, for(r=1,NQ+1, bigm[r,nv+j] = -mat3[r,j]));
    my(ker = matker(bigm));
    for(kcol=1,matsize(ker)[2],
      my(cv = vector(nv,i,ker[i,kcol]));
      if(cv == vector(nv), next);
      my(fq = O(x^(NQ+1))); for(i=1,nv, if(cv[i]!=0, fq += cv[i]*vec2ser(qexpv(mf1,V1[,i]~,NQ))));
      my(a0=polcoeff(fq,0)); if(a0==0, next); fq = fq/a0;
      my(qt = serreverse(tt));
      my(av = coefv(subst(fq,x,qt), NQ-6));
      my(intA=1); for(i=1,#av, if(denominator(av[i])!=1, intA=0; break));
      if(!intA, next);
      my(nl = #av-1);
      my(ne = matsize(V3e)[2]);
      my(bvs = vector(ne));
      for(i=1,ne,
        my(ph = vec2ser(qexpv(mf3e, V3e[,i]~, NQ)));
        bvs[i] = coefv(subst(fq*Dm2(ph - polcoeff(ph,0)), x, qt), nl);
      );
      my(good = select(m->(av[m+1]!=0), vector(nl,m,m)));
      my(n1 = good[#good\2], n3 = good[#good]);
      /* step 1: exact condition a_0(Phi) = 0 */
      my(a0v = vector(ne,i, polcoeff(vec2ser(qexpv(mf3e, V3e[,i]~, 4)),0)));
      my(k0 = matker(Mat(a0v)));
      my(ne0 = matsize(k0)[2]);
      /* step 2: within it, kill the slow direction (one numerical condition) */
      my(w1 = vector(ne0,j, sum(i=1,ne,k0[i,j]*bvs[i][n1+1])/av[n1+1]*1.0));
      my(w3 = vector(ne0,j, sum(i=1,ne,k0[i,j]*bvs[i][n3+1])/av[n3+1]*1.0));
      my(dv2 = w3 - w1);
      my(nrm = vecmax(vector(ne0,i,abs(dv2[i]))));
      my(k1 = if(nrm < 1e-30, matid(ne0), matker(Mat(vector(ne0,i, bestappr(dv2[i]/nrm, 10^10))))));
      my(dker = k0*k1);
      print("\n### ", tag, " pole-orb ",lb," lam1~", abs(av[n3+1]*1.)^(1/n3),
            "  dim Eis-eigsp=",ne, "  dim fold-regular=", matsize(dker)[2]);
      print("    A = ", vector(8,i,av[i]));
      for(j=1,matsize(dker)[2],
        my(c = dker[,j]);
        my(piv=0); for(i=1,ne, if(abs(c[i])>1e-20, piv=i));
        my(cq = vector(ne,i, bestappr(c[i]/c[piv],10^8)));
        my(bb = vector(nl+1, m, sum(i=1,ne, cq[i]*bvs[i][m])));
        my(xi = bb[n3+1]/av[n3+1]*1.0);
        my(dec = vector(3,jj, my(m=good[(jj*#good)\3]); my(z=bb[m+1]-xi*av[m+1]); if(z!=0, abs(z*1.)^(1/m),0)));
        my(kx=0); for(kz=1,4, my(ok=1); for(m=1,min(80,nl), if(denominator(dnn(m)^kz*bb[m+1])!=1, ok=0;break)); if(ok, kx=kz; break));
        print("    FR combo ", cq, "  xi = ", xi);
        print("        lindep(1,z2,G;xi) = ", lindep([1., Z2, GG, xi], 30), "   k=",kx, "   |b_n-xi a_n|^(1/n) = ", dec);
      );
    );
  ));
}
{RH = [ [12,[],[],"N=12"], [16,[16],[-I],"N=16+16"], [20,[4],[I],"N=20+4"],
        [36,[4,9],[I,1],"N=36+4,9"], [36,[4],[I],"N=36+4"] ];}
{for(h=1,#RH, analyse(RH[h][1],RH[h][2],RH[h][3],RH[h][4]));}
quit;
