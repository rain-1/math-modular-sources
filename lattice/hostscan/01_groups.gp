/* 01_groups.gp -- genus-0 groups Gamma_0(N)+W_S for N <= 120, with orbifold data.
   genus(X_0(N)/W_S) = dim of joint (+1)-eigenspace of {W_Q: Q in S} on S_2(Gamma_0(N)).
   Area: mu' = mu/|W_S|;  sum_i (1 - 1/m_i) = mu'/6 + 2 - nu_inf'   (genus 0).
   Reports: N, S, mu', nu_inf(Gamma_0(N)), nu_inf', A = sum(1-1/m), nu2, nu3 of Gamma_0(N).
   #Sigma = nu_inf' + #cone points.   Cone points all order 2  <=>  #cone = 2A. */
default(parisizemax, 6000000000);

idx(n) = my(f=factor(n)~); n*prod(i=1,#f, 1+1/f[1,i]);
nu2(n) = if(n%4==0, 0, my(f=factor(n)[,1]~); prod(i=1,#f, 1 + kronecker(-4, f[i])));
nu3(n) = if(n%9==0, 0, my(f=factor(n)[,1]~); prod(i=1,#f, 1 + kronecker(-3, f[i])));
ncusp0(n) = sumdiv(n, d, eulerphi(gcd(d,n/d)));

{almat(n,qq) = my(d=bezout(qq,n/qq)); [qq*d[1], 1; -n*d[2], qq];}

{cusplist(n) = my(res=List());
  fordiv(n,c, my(g=gcd(c,n/c)); my(seen=List());
    for(a=1,max(g,1),
      if(gcd(a,g)==1,
        my(new=1);
        for(k=1,#seen, if((a - seen[k])%g==0, new=0; break));
        if(new, listput(seen,a); listput(res,[a,c])))));
  Vec(res);}

cuspeq(n,p1,p2) = (p1[1]*p2[2] - p2[1]*p1[2]) % gcd(p1[2]*p2[2], n) == 0;

{cuspact(n,qq,p) = my(w=almat(n,qq));
  my(a = w[1,1]*p[1] + w[1,2]*p[2], c = w[2,1]*p[1] + w[2,2]*p[2]);
  my(g=gcd(a,c)); if(g!=0, a/=g; c/=g); if(c<0, a=-a; c=-c);
  if(c==0, [1,0], [a,c]);}

{norbits(n, invs) = my(cl=cusplist(n), nc=#cl, lab=vector(nc,i,i));
  for(i=1,nc, for(j=1,#invs,
    my(im = cuspact(n,invs[j],cl[i]));
    for(k=1,nc, if(cuspeq(n,im,cl[k]),
      my(o1=lab[i], o2=lab[k], m=min(o1,o2));
      for(z=1,nc, if(lab[z]==o1 || lab[z]==o2, lab[z]=m)); break))));
  #Set(Vec(lab));}

{allinv(qs) = my(r=#qs, res=List());
  for(mask=1, 2^r-1, my(q=1); for(j=1,r, if(bittest(mask,j-1), q*=qs[j])); listput(res,q));
  Vec(res);}

print("N|S|mu'|ncusp0|ncuspQ|A|nu2|nu3|Sigma_if_all2");
{for(nn=1,120,
  my(mf = mfinit([nn,2],1), d = mfdim(mf));
  my(gens = select(q->(q>1 && isprimepower(q) && gcd(q,nn/q)==1), divisors(nn)));
  my(mats = vector(#gens, i, if(d>0, my(al=mfatkininit(mf,gens[i])); al[2]/al[3], matrix(0,0))));
  for(mask=0, 2^#gens-1,
    my(sel = select(i->bittest(mask,i-1), vector(#gens,i,i)));
    my(g, qs = vector(#sel,j,gens[sel[j]]));
    if(d==0, g=0,
      my(P = matid(d));
      for(j=1,#sel, P = P*(matid(d)+mats[sel[j]])/2);
      g = matrank(P));
    if(g==0,
      my(ww = 2^#qs, mup = idx(nn)/ww);
      my(invs = if(#qs==0, [], allinv(qs)));
      my(nc0 = ncusp0(nn), ncq = norbits(nn, invs));
      my(A = mup/6 + 2 - ncq);
      print(nn,"|",qs,"|",mup,"|",nc0,"|",ncq,"|",A,"|",nu2(nn),"|",nu3(nn),"|",ncq+2*A);
    )));}
quit;
