/* 06_rows.gp -- build the chi_{-4} rows on every genus-0 host Gamma_0(N)+W
   (4|N) whose Atkin-Lehner eigenspaces are rational, i.e. every Q in W is a
   perfect square.  For each choice of the cusp orbit carrying the pole of the
   hauptmodul t (t = q + O(q^2), simple zero at infinity), we
     - build t as a (root of a) modular unit via Ligozat,
     - solve linearly for the weight-1 form F with F * D t in M_3(N,chi_{-4}),
     - form A(t) = F o q(t) and test integrality,
     - fit the minimal recurrence and read off the singular set.
*/
default(parisizemax, 8000000000);
default(realprecision, 60);
read("lattice/catalan_al_hosts/lib2.gp");

NQ = 90;

{RH = [ [4,[],[]], [4,[4],[-I]], [8,[],[]], [12,[],[]], [12,[4],[I]], [12,[4],[-I]],
        [16,[],[]], [16,[16],[I]], [16,[16],[-I]], [20,[4],[I]],
        [36,[4],[I]], [36,[4],[-I]], [36,[4,9],[I,1]], [36,[4,9],[I,-1]], [36,[4,9],[-I,-1]] ];}

aldiv(n,qq,d) = (qq/gcd(d,qq))*gcd(d,n/qq);
{allinv(qs) = my(r=#qs, res=List(1)); for(mask=1, 2^r-1, my(q=1); for(j=1,r, if(bittest(mask,j-1), q*=qs[j])); listput(res,q)); Vec(res);}
ncu(n,c) = eulerphi(gcd(c,n/c));

{ratproj(mm, ev) =
  my(d=matsize(mm)[1]);
  my(others = select(z->abs(z-ev)>1e-10, [1,-1,I,-I]));
  my(p1 = matid(d)); for(k=1,#others, p1 = p1*(mm - others[k]*matid(d))/(ev-others[k]));
  matrix(d,d,i,j, bestappr(real(p1[i,j]),10^6));
}

{qexpv(mfx, cv, nn) = my(bb=mfbasis(mfx), s=vector(nn+1));
  for(i=1,#bb, if(cv[i]!=0, s += cv[i]*mfcoefs(bb[i],nn))); s;}
vec2ser(v) = sum(n=1,#v, v[n]*x^(n-1)) + O(x^#v);
coefv(s,nn) = vector(nn+1, i, polcoeff(s,i-1));
Dop(s) = my(n=serprec(s,x)); sum(k=1,n-1, k*polcoeff(s,k)*x^k) + O(x^n);

{run(nn, qs, evs) =
  my(dv = divisors(nn), kk = #dv);
  my(invs = allinv(qs));
  /* orbits of W on cusp-denominators */
  my(orb = vector(kk,i,i));
  for(i=1,kk, for(j=1,#invs, my(d2 = aldiv(nn,invs[j],dv[i]));
    for(l=1,kk, if(dv[l]==d2, my(m=min(orb[i],orb[l]), o1=orb[i], o2=orb[l]);
      for(z=1,kk, if(orb[z]==o1||orb[z]==o2, orb[z]=m))))));
  my(oinf = orb[kk]);   /* orbit label of the cusp of denominator N = infinity */
  my(labs = Set(Vec(orb)));
  my(mf1 = mfinit([nn,1,-4],4), mf3 = mfinit([nn,3,-4],4));
  my(d1 = mfdim(mf1), d3 = mfdim(mf3));
  my(p1 = matid(d1), p3 = matid(d3));
  for(j=1,#qs, my(al1=mfatkininit(mf1,qs[j]), al3=mfatkininit(mf3,qs[j]));
    p1 = p1*ratproj(al1[2]/al1[3], evs[j]); p3 = p3*ratproj(al3[2]/al3[3], evs[j]));
  my(V1 = matimage(p1), V3 = matimage(p3));
  print("\n===== N=",nn," W=",qs," ev=",evs,"  dimV1=",matsize(V1)[2]," dimV3=",matsize(V3)[2]);
  print("      cusp-denominator orbits: ", vector(kk,i,[dv[i],orb[i]]));
  my(b3 = mfbasis(mf3));
  my(mat3 = matrix(NQ+1, d3, i,j, mfcoefs(b3[j],NQ)[i]));
  my(b1 = mfbasis(mf1));
  for(li=1,#labs, my(lb=labs[li]); if(lb != oinf,
    /* ramification index of a cusp = |W| / (orbit size in cusps) */
    my(ninf = sum(i=1,kk, if(orb[i]==oinf, ncu(nn,dv[i]), 0)));
    my(npol = sum(i=1,kk, if(orb[i]==lb,   ncu(nn,dv[i]), 0)));
    my(einf = #invs/ninf, epol = #invs/npol);
    if(denominator(einf)!=1 || denominator(epol)!=1, next);
    my(ordv = vector(kk));
    for(i=1,kk, if(orb[i]==oinf, ordv[i]=einf); if(orb[i]==lb, ordv[i]=-epol));
    my(rr = unitexp(nn, ordv));
    my(sr = sum(i=1,kk, rr[i]));
    if(sr != 0, print("   pole-orbit ",lb,": sum r = ",sr," (not weight 0) -- skip"); next);
    my(tt = hauptmod(nn, rr, NQ+1));
    /* normalise t = q + O(q^2) : hauptmod already gives q*exp(...) */
    my(dt = Dop(tt));
    /* solve for F in V1 with F*Dt in M_3 */
    my(nv = matsize(V1)[2], rows = List());
    my(cand = List());
    for(i=1,nv,
      my(fq = vec2ser(qexpv(mf1, V1[,i]~, NQ)));
      listput(cand, coefv(fq*dt + O(x^(NQ+1)), NQ));
    );
    /* find combos c with sum c_i cand_i in column space of mat3 */
    my(bigm = matrix(NQ+1, nv+d3));
    for(i=1,nv, for(r=1,NQ+1, bigm[r,i] = cand[i][r]));
    for(j=1,d3, for(r=1,NQ+1, bigm[r,nv+j] = -mat3[r,j]));
    my(ker = matker(bigm));
    my(good = List());
    for(kcol=1,matsize(ker)[2],
      my(cv = vector(nv,i,ker[i,kcol]));
      if(cv != vector(nv), listput(good, cv)));
    print("   pole-orbit ",lb,":  t = ", tt + O(x^8), "   #F-sol=", #good);
    if(#good==0, next);
    
    for(g=1,#good,
      my(cv = good[g]);
      my(fq = O(x^(NQ+1))); for(i=1,nv, if(cv[i]!=0, fq += cv[i]*vec2ser(qexpv(mf1,V1[,i]~,NQ))));
      my(a0 = polcoeff(fq,0)); if(a0==0, print("      F has a_0=0 -- skip"); next);
      fq = fq/a0;
      my(qt = serreverse(tt));
      my(av = coefv(subst(fq,x,qt), min(60,NQ-2)));
      my(intA = 1); for(i=1,min(#av,50), if(denominator(av[i])!=1, intA=0; break));
      print("      F = ", fq+O(x^7));
      print("      A = ", vector(min(9,#av),i,av[i]), if(intA," [integral]"," [NOT integral]"));
    );
  ));
}

{for(h=1,#RH, run(RH[h][1], RH[h][2], RH[h][3]));}
quit;
