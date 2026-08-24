/* Gamma_0(12) chi_{-4} rows, pole at cusp c=3 ("placement 3") and c=6 ("placement 6").
   Dumps A(t) and the Eisenstein companions B_Phi(t). */
default(parisize, 8000000000);
default(realprecision, 120);
read("/home/ubuntu/code/math-modular-sources/lattice/catalan_al_hosts/lib2.gp");
NQ = 300;
OUT = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/";
vec2ser(v) = sum(n=1,#v, v[n]*x^(n-1)) + O(x^#v);
coefv(s,nn) = vector(nn+1, i, polcoeff(s,i-1));
Dop(s) = my(n=serprec(s,x)); sum(k=1,n-1, k*polcoeff(s,k)*x^k) + O(x^n);
Dm2(s) = my(n=serprec(s,x)); sum(k=1,n-1, polcoeff(s,k)/k^2*x^k) + O(x^n);
dn(n) = lcm(vector(max(n,1),j,j));
NN = 12; DV = divisors(NN);   /* [1,2,3,4,6,12] */
mf1 = mfinit([NN,1,-4],4); mf3 = mfinit([NN,3,-4],4); mf3e = mfinit([NN,3,-4],3);
d1 = mfdim(mf1); d3 = mfdim(mf3); d3e = mfdim(mf3e);
print("dims: M1=",d1," M3=",d3," Eis3=",d3e);
b1 = mfbasis(mf1); b3 = mfbasis(mf3); b3e = mfbasis(mf3e);
mat3 = matrix(NQ+1, d3, i,j, mfcoefs(b3[j],NQ)[i]);
{doplace(cidx, tag) =
  my(ordv = vector(#DV));
  ordv[#DV] = 1; ordv[cidx] = -1;
  my(rr = unitexp(NN, ordv));
  print("\n===== ", tag, "  pole at cusp c=", DV[cidx], "   eta-exponents r = ", rr, "  sum=", sum(i=1,#rr,rr[i]));
  my(tt = hauptmod(NN, rr, NQ+1), dt = Dop(tt));
  print("  t = ", tt + O(x^8));
  my(cand = List());
  for(i=1,d1, listput(cand, coefv(vec2ser(mfcoefs(b1[i],NQ))*dt + O(x^(NQ+1)), NQ)));
  my(bigm = matrix(NQ+1, d1+d3));
  for(i=1,d1, for(r=1,NQ+1, bigm[r,i] = cand[i][r]));
  for(j=1,d3, for(r=1,NQ+1, bigm[r,d1+j] = -mat3[r,j]));
  my(ker = matker(bigm));
  print("  #F solutions: ", matsize(ker)[2]);
  for(kcol=1,matsize(ker)[2],
    my(cv = vector(d1,i,ker[i,kcol]));
    if(cv == vector(d1), next);
    my(fq = O(x^(NQ+1))); for(i=1,d1, if(cv[i]!=0, fq += cv[i]*vec2ser(mfcoefs(b1[i],NQ))));
    my(a0=polcoeff(fq,0)); if(a0==0, next); fq = fq/a0;
    my(qt = serreverse(tt), NM = NQ-4);
    my(av = coefv(subst(fq,x,qt), NM));
    my(intA=1, badn=0); for(i=1,#av, if(denominator(av[i])!=1, intA=0; badn=i-1; break));
    print("  F = ", fq + O(x^7));
    print("  A = ", vector(12,i,av[i]));
    print("  a_n integral to n=",NM,"? ", intA, if(badn, Str(" first non-integer n=",badn),""));
    my(fn = Str(OUT,"lvl12_p",DV[cidx],"_A.txt")); write1(fn,"");
    for(i=1,#av, write(fn, av[i]));
    print("  wrote ", fn);
    /* Eisenstein companions */
    for(i=1,d3e,
      my(ph = vec2ser(mfcoefs(b3e[i],NQ)), a00 = polcoeff(ph,0));
      my(bs = subst(fq*Dm2(ph - a00), x, qt));
      my(bv = coefv(bs, NM));
      my(kk2=0); for(kx=0,4, my(ok=1); for(m=1,min(60,#bv-1), if(denominator(dn(m)^kx*bv[m+1])!=1, ok=0;break)); if(ok, kk2=kx; break));
      my(rat = bv[NM+1]/av[NM+1]*1.0);
      print("   Phi_",i," (a_0=",a00,")  b_n/a_n = ", rat, "   k=", kk2);
      print("      lindep([1,zeta(2),G,ratio]) = ", lindep([1., zeta(2), Catalan, rat], 40));
      my(fnb = Str(OUT,"lvl12_p",DV[cidx],"_B",i,".txt")); write1(fnb,"");
      for(m=1,#bv, write(fnb, bv[m]));
    );
  );
}
doplace(3, "placement 3");
doplace(5, "placement 6");
quit;
