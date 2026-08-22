/* 20_level8_slopes.gp -- measured 2-adic slopes in the symmetrised coordinate
   on the LEVEL-8 Catalan host (Zagier E).

   Host: t = q eta_1^4 eta_4^2 eta_8^4 / eta_2^10,  F = eta_2^10/(eta_1^4 eta_4^4),
   source Phi_E = (1-8V2)E,  lambda_1=8 (fold t=1/8), lambda_2=4, s=1/lambda_2=1/4.
   CDT descent: sigma(t) = s t/(t-s) = t/(4t-1),  y = t + sigma(t) = 4t^2/(4t-1).

   Because s>0 here, t(v) with y=v^2 is not real; we therefore extract the
   y-expansion of a sigma-invariant series directly in the t-variable
   (the `to_y` algorithm of lattice/cdt_finder/indep_check2.py): y has t-order 2,
   so the coefficients are read off orders t^0, t^2, t^4, ... successively.
   This is exact over Q and works for either sign of s.                        */
default(parisizemax, 12*10^9);
NTERM = 302;
read("/home/ubuntu/code/math-modular-sources/lattice/catalan_two_classes/00_setup.gp");
H8 = mkhost(t8,F8); Av = Aof(H8); NMAX = NT-2;
BE8 = Bof(H8, mkPhi(Sin,[[1,1],[2,-8]]));      /* canonical (1-8V2)E, xi=G/2   */
B0  = Bof(H8, mkPhi(Sin,[[1,1],[2,-4]]));      /* target zero (1-4V2)E          */
sig(f) = subst(f, q, q/(4*q-1) + O(q^(NMAX+1)));
yx = 4*q^2/(4*q-1) + O(q^(NMAX+1));            /* y(t), t-order 2               */
KY = NMAX\2 - 3;
{toy(gx) =                                     /* gx sigma-invariant -> series in y */
  my(cur = gx, r = 0, yp = 1 + O(q^(NMAX+1)));
  for(n=0,KY,
    my(c = polcoeff(cur,2*n)/polcoeff(yp,2*n));
    r += c*q^n; cur -= c*yp; yp *= yx);
  r + O(q^(KY+1));}
{toyodd(gx) =                                  /* gx sigma-anti-invariant: divide by t-sigma(t) */
  my(od = q - sig(q));
  toy(gx/od);}
rate(sr,n) = abs(polcoeff(sr,n))^(1./n);
v2s(c,k) = if(c==0, "zero", valuation(c,2)/(k*1.));
{rp(nm, f, odd) =
  my(g = if(odd, f - sig(f), f + sig(f)), Y = if(odd, toyodd(g), toy(g)));
  print("  ", nm, "  v2(d_k)/k at k=30,60,",KY,": ",
        [v2s(polcoeff(Y,30),30), v2s(polcoeff(Y,60),60), v2s(polcoeff(Y,KY),KY)],
        "   |d_k|^1/k: ", [rate(Y,60), rate(Y,KY)]);}
Li1 = sum(n=1,NMAX,4^n/n*q^n)+O(q^(NMAX+1));
Li2 = sum(n=1,NMAX,4^n/n^2*q^n)+O(q^(NMAX+1));
Li3 = sum(n=1,NMAX,4^n/n^3*q^n)+O(q^(NMAX+1));
print("LEVEL-8 host: t = ", t8 + O(q^7));
print("A(t) = ", Av + O(q^7));
print("B_E(t) = ", BE8 + O(q^7));
print("sanity: sigma(sigma(t)) - t = ", sig(sig(q)) - q + O(q^6));
print("sanity: y(t) - (t + sigma(t)) = ", yx - (q + sig(q)) + O(q^8));
print("\n=== 2-adic slopes in y = 4t^2/(4t-1),  k <= ", KY, " (t-series to n=", NMAX, ") ===");
rp("Sym+ 1                      ", 1+O(q^(NMAX+1)), 0);
rp("Sym+ Li1(4t)                ", Li1, 0);
rp("Sym+ Li2(4t)                ", Li2, 0);
rp("Sym+ Li3(4t)                ", Li3, 0);
rp("Sym- Li2(4t) / (t-sigma t)  ", Li2, 1);
rp("Sym- Li3(4t) / (t-sigma t)  ", Li3, 1);
rp("Sym+ A            (host)    ", Av, 0);
rp("Sym- A / (t-sigma t)        ", Av, 1);
rp("Sym+ B_E          (raw)     ", BE8, 0);
rp("Sym+ B_E + (7/3)A (cond)    ", BE8 + 7/3*Av, 0);
rp("Sym+ B_E + (1/8)A (cond)    ", BE8 + 1/8*Av, 0);
rp("Sym+ B_E - (1/2)A (cond)    ", BE8 - 1/2*Av, 0);
rp("Sym- (B_E+(7/3)A)/(t-sig t) ", BE8 + 7/3*Av, 1);
rp("Sym+ B_0 = (1-4V2)E cmpn    ", B0, 0);
print("\n(for reference, slopes in the t-coordinate)");
{for(j=1,4, my(f=[Av,BE8,B0,Li2][j], nm=["A       ","B_E     ","B_0     ","Li2(4t) "][j]);
  print("  ", nm, " v2(c_n)/n at n=100,200,",NMAX,": ",
     [v2s(polcoeff(f,100),100), v2s(polcoeff(f,200),200), v2s(polcoeff(f,NMAX),NMAX)]));}
quit
