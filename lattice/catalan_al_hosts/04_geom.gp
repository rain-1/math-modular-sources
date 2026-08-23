/* 04_geom.gp -- orbifold data of the genus-0 hosts Gamma_0(N)+W, 4 | N.
   Since 4 | N one has nu_2(Gamma_0(N)) = nu_3(Gamma_0(N)) = 0, so every
   elliptic point of the quotient is of order 2 and comes from an interior
   fixed point of an Atkin-Lehner involution.
   Cusp action: W_Q = [Q*aa, bb; N*cc, Q*dd] with Q*aa*dd - (N/Q)*bb*cc = 1.
   Cusp equivalence (Diamond-Shurman 3.8.3): a1/c1 ~ a2/c2 mod Gamma_0(N)
   iff a1*c2 = a2*c1 (mod gcd(c1*c2, N)).
*/
default(parisizemax, 8000000000);

{HOSTS = [ [4,[]], [4,[4]], [8,[]], [8,[8]], [12,[]], [12,[3]], [12,[4]], [12,[3,4]],
          [16,[]], [16,[16]], [20,[4]], [20,[4,5]], [24,[8]], [24,[3,8]],
          [28,[7]], [28,[4,7]], [32,[32]], [36,[4]], [36,[4,9]],
          [44,[4,11]], [56,[7,8]], [60,[3,4,5]], [92,[4,23]] ];}

idx(n) = my(f=factor(n)~); n*prod(i=1,#f, 1+1/f[1,i]);

{almat(n,qq) = my(g,u,v); my(d=bezout(qq,n/qq)); /* d=[u,v,g] with u*qq+v*(n/qq)=1 */
  my(aa=d[1], bb=1, cc=-d[2], dd=1);
  /* want Q*aa*dd - (N/Q)*bb*cc = 1 : Q*aa + (N/Q)*d[2] = 1 -> set bb=1, cc=-d[2] */
  [qq*aa, bb; n*cc, qq*dd];}

/* list of cusps of Gamma_0(N) as [a,c] representatives */
{cusplist(n) = my(res=List());
  fordiv(n,c, my(g=gcd(c,n/c));
    my(seen=List());
    for(a=1,max(g,1), if(gcd(a,c)==1 || c==1,
      if(gcd(a,g)==1,
        my(new=1);
        for(k=1,#seen, if((a*c - seen[k]*c)%gcd(c*c,n)==0, new=0; break));
        if(new, listput(seen,a); listput(res,[a,c]))))));
  Vec(res);}

cuspeq(n,p1,p2) = (p1[1]*p2[2] - p2[1]*p1[2]) % gcd(p1[2]*p2[2], n) == 0;

{cuspact(n,qq,p) = my(w=almat(n,qq));
  my(a = w[1,1]*p[1] + w[1,2]*p[2], c = w[2,1]*p[1] + w[2,2]*p[2]);
  my(g=gcd(a,c)); if(g!=0, a/=g; c/=g); if(c<0, a=-a; c=-c);
  if(c==0, [1,0], my(cc=gcd(c,n)); [a,c]);}

{fixpts(n, qq) = my(mf=mfinit([n,2],1), d=mfdim(mf)); if(d==0, return(2));
  my(al=mfatkininit(mf,qq)); round(2 - 2*trace(al[2]/al[3]));}
{allinv(qs) = my(r=#qs, res=List());
  for(mask=1, 2^r-1, my(q=1); for(j=1,r, if(bittest(mask,j-1), q*=qs[j])); listput(res,q));
  Vec(res);}

print("N  W          mu'    ncusp0 ncuspG  e2  #sing   (area check)");
{for(h=1,#HOSTS,
  my(nn=HOSTS[h][1], qs=HOSTS[h][2], ww=2^#qs);
  my(mu=idx(nn), mup=mu/ww);
  my(cl = cusplist(nn), nc0 = #cl);
  my(invs = if(#qs==0, [], allinv(qs)));
  /* orbits of W on cusps */
  my(lab = vector(nc0,i,i));
  for(i=1,nc0, for(j=1,#invs,
    my(im = cuspact(nn,invs[j],cl[i]));
    for(k=1,nc0, if(cuspeq(nn,im,cl[k]), my(m=min(lab[i],lab[k])); my(o1=lab[i],o2=lab[k]);
      for(z=1,nc0, if(lab[z]==o1 || lab[z]==o2, lab[z]=m)); break))));
  my(norb = #Set(Vec(lab)));
  my(e2 = 4 + mup/3 - 2*norb);
  print(nn," ",qs,"\t",mup,"\t",nc0,"\t",norb,"\t",e2,"\t",norb+e2);
);}
quit;
