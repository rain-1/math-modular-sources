/* 07_sing.gp -- for every row found by 06_rows.gp, fit the second-order Fuchsian
   ODE  q2(t) A'' + q1(t) A' + q0(t) A = 0  and read the singular set off q2.
   Also reports lambda_1 = 1/|t_1| (fold = nearest singular point to 0) and
   lambda_2 = 1/|t_2|. */
default(parisizemax, 8000000000);
default(realprecision, 60);
read("lattice/catalan_al_hosts/lib2.gp");

NQ = 200;
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

/* fit q2 A'' + q1 A' + q0 A = 0, deg q_i <= dg, from the coefficient vector av */
{fitode(av, dg) =
  my(nn = #av-1);
  my(as = sum(i=0,nn, av[i+1]*x^i) + O(x^(nn+1)));
  my(a1 = deriv(as), a2 = deriv(a1));
  my(use = nn-4);
  my(cols = List());
  for(e=0,dg, listput(cols, coefv(x^e*a2 + O(x^(use+1)), use)));
  for(e=0,dg, listput(cols, coefv(x^e*a1 + O(x^(use+1)), use)));
  for(e=0,dg, listput(cols, coefv(x^e*as + O(x^(use+1)), use)));
  my(nc = #cols, mm = matrix(use+1, nc, i,j, cols[j][i]));
  my(ker = matker(mm));
  if(matsize(ker)[2]==0, return(0));
  my(kv = ker[,1]);
  sum(e=0,dg, kv[e+1]*x^e);
}

{report(nn,qs,evs,tt,fq,tag) =
  my(qt = serreverse(tt));
  my(av = coefv(subst(fq,x,qt), NQ-4));
  my(intA=1); for(i=1,#av, if(denominator(av[i])!=1, intA=0; break));
  my(q2=0, dg=0);
  for(d=1,10, q2 = fitode(av, d); if(q2!=0 && poldegree(q2)>0, dg=d; break));
  if(q2==0, print("   ",tag,": no ODE of degree <= 10"); return(0));
  q2 = q2/content(q2);
  my(rts = polroots(q2));
  my(fin = select(z->abs(z)>1e-25, rts));
  my(mods = vecsort(vector(#fin,i,abs(fin[i]))));
  print("   ",tag," intA=",intA," degq2=",poldegree(q2), "  q2 = ", q2);
  print("        roots: ", vector(#rts,i, [bestappr(real(rts[i]),10^8), bestappr(imag(rts[i]),10^8)]));
  print("        |t| sorted: ", vector(#mods,i,mods[i]));
  if(#mods>=1, print("        lambda_1 = ", 1/mods[1], if(#mods>=2, Str("   lambda_2 = ", 1/mods[2]), "")));
  [av,q2];
}

{run(nn, qs, evs) =
  my(dv = divisors(nn), kk = #dv, invs = allinv(qs));
  my(orb = vector(kk,i,i));
  for(i=1,kk, for(j=1,#invs, my(d2 = aldiv(nn,invs[j],dv[i]));
    for(l=1,kk, if(dv[l]==d2, my(m=min(orb[i],orb[l]), o1=orb[i], o2=orb[l]);
      for(z=1,kk, if(orb[z]==o1||orb[z]==o2, orb[z]=m))))));
  my(oinf = orb[kk], labs = Set(Vec(orb)));
  my(mf1 = mfinit([nn,1,-4],4), mf3 = mfinit([nn,3,-4],4));
  my(d1 = mfdim(mf1), d3 = mfdim(mf3));
  my(p1 = matid(d1), p3 = matid(d3));
  for(j=1,#qs, my(al1=mfatkininit(mf1,qs[j]), al3=mfatkininit(mf3,qs[j]));
    p1 = p1*ratproj(al1[2]/al1[3], evs[j]); p3 = p3*ratproj(al3[2]/al3[3], evs[j]));
  my(V1 = matimage(p1), V3 = matimage(p3));
  print("\n===== N=",nn," W=",qs," ev=",evs);
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
    my(tt = hauptmod(nn, rr, NQ+1), dt);  dt = Dop(tt);
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
      report(nn,qs,evs,tt,fq, Str("pole-orb ",lb," F=",fq+O(x^5)));
    );
  ));
}

{RH = [ [4,[],[]], [4,[4],[-I]], [8,[],[]], [12,[],[]],
        [16,[],[]], [16,[16],[I]], [16,[16],[-I]], [20,[4],[I]],
        [36,[4],[I]], [36,[4],[-I]], [36,[4,9],[I,1]], [36,[4,9],[I,-1]], [36,[4,9],[-I,-1]] ];}
{for(h=1,#RH, run(RH[h][1], RH[h][2], RH[h][3]));}
quit;
