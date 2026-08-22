/* Level-16 host: (i) exact fold-regularity functional 8c2+c4=0;
   (ii) the normaliser descent sigma(x) = -x/(4x+1), y = 4x^2/(4x+1);
   (iii) 2-adic slopes of the pure module and of the doubly-small function in y. */
read("00_setup.gp");
G = Catalan;
H16 = mkhost(x16,F16); Av = Aof(H16); NMAX = NT-2;
BE1 = Bof(H16, mkPhi(Sin,[[1,1]]));
BE2 = Bof(H16, mkPhi(Sin,[[2,1]]));
BE4 = Bof(H16, mkPhi(Sin,[[4,1]]));
rate(S,n) = abs(polcoeff(S,n))^(1./n);
P2(c1,c2,c4) = c1+c2/4+c4/16;
{ TESTS = [[1,1,-8],[1,-1,8],[2,3,-24],[1,0,0],[0,1,-8],[1,4,-32],[1,1,1],[1,0,-8]]; }
print("=== fold-regularity test: ell = 8*c2+c4 ===");
{for(i=1,#TESTS, my(c=TESTS[i], Bv=c[1]*BE1+c[2]*BE2+c[3]*BE4, xi=-P2(c[1],c[2],c[3])*G/2, Cv=Bv-xi*Av);
  print("  (c1,c2,c4)=",c,"  ell=",8*c[2]+c[3],
        "  b/a(n=",NMAX,")=",polcoeff(Bv,NMAX)/polcoeff(Av,NMAX)*1.," pred ",xi,
        "  |c_n|^1/n=",rate(Cv,NMAX)));}

/* ---------- normaliser descent ---------- */
sig(f) = subst(f, q, -q/(4*q+1) + O(q^(NMAX+1)));
v = 'v;
xofv = (v^2 + v*sqrt(1+v^2+O(v^(2*NMAX))))/2;          /* x(v), y=v^2 */
toY(g) = my(S=subst(g,q,xofv)); S;   /* series in v; should be even */
{evencheck(S,nm) = my(bad=-1); for(k=1,NMAX, if(k%2==1 && polcoeff(S,k)!=0, bad=k;break)); print("   ",nm," first odd nonzero v-coeff: ",bad);}
inY(S) = my(r=0); for(k=0,NMAX\2, r += polcoeff(S,2*k)*q^k); r + O(q^(NMAX\2+1));

B0 = BE1+4*BE2-32*BE4;
print("\n=== 2-adic slopes in x ===");
print("  B0 : v2/n at n=60,100,",NMAX,": ",[valuation(polcoeff(B0,60),2)/60.,valuation(polcoeff(B0,100),2)/100.,valuation(polcoeff(B0,NMAX),2)/(NMAX*1.)]);
Li2 = sum(n=1,NMAX,(-4)^n/n^2*q^n)+O(q^(NMAX+1));
Li1 = sum(n=1,NMAX,(-4)^n/n*q^n)+O(q^(NMAX+1));
print("  Li2(-4x): v2/n at n=60,100: ",[valuation(polcoeff(Li2,60),2)/60.,valuation(polcoeff(Li2,100),2)/100.]);

print("\n=== descend to y ===");
{ for(j=1,4,
   my(f = [B0, Li2, Li1, Av][j], nm = ["B0","Li2(-4x)","Li1(-4x)","A"][j]);
   my(sp = f + sig(f), sm = f - sig(f));
   my(SP = toY(sp), SM = toY(sm));
   evencheck(SP, concat("Sym+ ",nm));
   my(YP = inY(SP));
   my(K = NMAX\2 - 1);
   print("   Sym+ ",nm,": in y, v2(d_k)/k at k=20,40,",K,": ",
     [valuation(polcoeff(YP,20),2)/20., valuation(polcoeff(YP,40),2)/40., valuation(polcoeff(YP,K),2)/(K*1.)]);
   print("        |d_k|^{1/k} at k=20,40,",K,": ",[rate(YP,20),rate(YP,40),rate(YP,K)]);
 ); }
