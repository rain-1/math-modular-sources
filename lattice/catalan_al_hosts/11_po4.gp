/* 11_po4.gp -- the Gamma_0(12) hauptmodul placement with singular set
   {0, +-1/3, +-1, infty}: lambda_1 = 3 attained at BOTH t = 1/3 and t = -1/3,
   lambda_2 = 1.  Question: is there a rational combination of A and the four
   Eisenstein companions whose coefficients grow like 1^n (i.e. regular at both
   folds)?  Test by LLL on the coefficient vectors at large n. */
default(parisizemax, 8000000000);
default(realprecision, 60);
read("lattice/catalan_al_hosts/lib2.gp");
NQ = 220;
{qexpv(mfx, cv, nn) = my(bb=mfbasis(mfx), s=vector(nn+1));
  for(i=1,#bb, if(cv[i]!=0, s += cv[i]*mfcoefs(bb[i],nn))); s;}
vec2ser(v) = sum(n=1,#v, v[n]*x^(n-1)) + O(x^#v);
coefv(s,nn) = vector(nn+1, i, polcoeff(s,i-1));
Dop(s) = my(n=serprec(s,x)); sum(k=1,n-1, k*polcoeff(s,k)*x^k) + O(x^n);
Dm2(s) = my(n=serprec(s,x)); sum(k=1,n-1, polcoeff(s,k)/k^2*x^k) + O(x^n);

nn=12; dv=divisors(nn);
/* pole-orbit 4 = the cusp of denominator 4 */
ordv = vector(#dv); for(i=1,#dv, if(dv[i]==12, ordv[i]=1); if(dv[i]==4, ordv[i]=-1));
rr = unitexp(nn, ordv); print("eta exponents r = ", rr, "  sum=", sum(i=1,#dv,rr[i]));
tt = hauptmod(nn, rr, NQ+1); print("t = ", tt+O(x^9));
mf1 = mfinit([nn,1,-4],4); mf3 = mfinit([nn,3,-4],4); mf3e = mfinit([nn,3,-4],3);
dt = Dop(tt);
b1 = mfbasis(mf1); d1=#b1; b3=mfbasis(mf3); d3=#b3;
mat3 = matrix(NQ+1,d3,i,j,mfcoefs(b3[j],NQ)[i]);
cand = vector(d1, i, coefv(vec2ser(mfcoefs(b1[i],NQ))*dt + O(x^(NQ+1)), NQ));
bigm = matrix(NQ+1, d1+d3);
for(i=1,d1, for(r=1,NQ+1, bigm[r,i]=cand[i][r]));
for(j=1,d3, for(r=1,NQ+1, bigm[r,d1+j]=-mat3[r,j]));
ker = matker(bigm);
print("F-solutions: ", matsize(ker)[2]);
{for(kc=1,matsize(ker)[2],
  my(cv=vector(d1,i,ker[i,kc])); if(cv==vector(d1), next);
  my(fq=O(x^(NQ+1))); for(i=1,d1, if(cv[i]!=0, fq += cv[i]*vec2ser(mfcoefs(b1[i],NQ))));
  my(a0=polcoeff(fq,0)); if(a0==0, next); fq=fq/a0;
  my(qt=serreverse(tt));
  my(av=coefv(subst(fq,x,qt), NQ-6));
  print("\nF = ", fq+O(x^7)); print("A = ", vector(9,i,av[i]));
  my(b3e=mfbasis(mf3e), ne=#b3e);
  my(bvs=vector(ne, i, my(ph=vec2ser(mfcoefs(b3e[i],NQ))); coefv(subst(fq*Dm2(ph-polcoeff(ph,0)),x,qt), NQ-6)));
  /* LLL: look for c_0 A + sum c_i B_i with coefficients of size O(1) */
  my(nl=#av-1);
  my(W = 12);                       /* window of consecutive n */
  my(n0 = nl-W+1);
  my(sc = 10^40);
  my(rows = ne+1);
  my(mm = matrix(rows, rows+W));
  for(i=1,rows, mm[i,i]=1);
  for(j=1,W, my(n=n0+j-1);
    my(den = abs(av[n+1]*1.));
    mm[1,rows+j] = round(sc*av[n+1]/den);
    for(i=1,ne, mm[i+1,rows+j] = round(sc*bvs[i][n+1]/den)));
  my(T = qflll(mm~));
  print("   LLL shortest combos (c_0=A, c_1..c_",ne,"):");
  for(k=1,min(3,matsize(T)[2]),
    my(c = T[,k]);
    my(cc = vector(rows,i,c[i]));
    my(gg = vector(4,jj, my(n=[nl\3, nl\2, (2*nl)\3, nl][jj]);
        my(z = cc[1]*av[n+1] + sum(i=1,ne, cc[i+1]*bvs[i][n+1]));
        if(z!=0, abs(z*1.)^(1/n), 0)));
    print("     c=",cc,"   |combo_n|^(1/n) = ", gg));
);}
quit;
