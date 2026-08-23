/* 08_full.gp -- full row data on every chi_{-4} host: A to high order, growth
   rate |a_n|^{1/n} (which singular points are genuine), the Eisenstein
   companions B_Phi = F D^{-2}Phi in t, their periods b_n/a_n identified against
   {1, zeta(2), G}, the sharp denominator exponent k, and the decay
   |b_n - xi a_n|^{1/n} of the fold-regular combinations. */
default(parisizemax, 8000000000);
default(realprecision, 60);
read("lattice/catalan_al_hosts/lib2.gp");

NQ = 240;
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
  my(b3e = mfbasis(mf3e));
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
      my(as = subst(fq,x,qt));
      my(av = coefv(as, NQ-4));
      my(intA=1); for(i=1,#av, if(denominator(av[i])!=1, intA=0; break));
      if(!intA, next);
      my(nl = #av-1);
      print("\n### ", tag, "  pole-orb ",lb,"  F=",fq+O(x^5));
      print("   A = ", vector(9,i,av[i]));
      print("   |a_n|^(1/n) at n=", [nl\4, nl\2, nl], ": ",
            vector(3,j, my(m=[nl\4,nl\2,nl][j]); if(av[m+1]!=0, abs(av[m+1]*1.)^(1/m), 0)));
      /* Eisenstein companions in the matching eigenspace */
      my(ne = matsize(V3e)[2]);
      print("   Eisenstein classes in the eigenspace: dim ", ne);
      for(i=1,ne,
        my(ph = vec2ser(qexpv(mf3e, V3e[,i]~, NQ)));
        my(a00 = polcoeff(ph,0));
        my(php = ph - a00);
        my(bs = subst(fq*Dm2(php), x, qt));
        my(bv = coefv(bs, NQ-4));
        my(kk2=0); for(kx=1,4, my(ok=1); for(m=1,min(60,#bv-1), if(denominator(dnn(m)^kx*bv[m+1])!=1, ok=0;break)); if(ok, kk2=kx; break));
        my(rat = if(av[nl+1]!=0, bv[nl+1]/av[nl+1]*1., 0));
        my(ld = lindep([1., Z2, GG, rat], 25));
        print("     Phi_",i," a_0=",a00,"  b_n/a_n(n=",nl,") = ", rat, "  k=",kk2, "  lindep(1,z2,G,r)=", ld);
      );
    );
  ));
}

{RH = [ [8,[],[],"N=8"], [12,[],[],"N=12"], [16,[],[],"N=16"],
        [16,[16],[-I],"N=16+16"], [20,[4],[I],"N=20+4"],
        [36,[4,9],[I,1],"N=36+4,9"] ];}
{for(h=1,#RH, analyse(RH[h][1],RH[h][2],RH[h][3],RH[h][4]));}
quit;
