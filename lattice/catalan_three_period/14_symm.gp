/* 14_symm.gp -- descend the level-16 objects along the normaliser involution
   sigma(x) = -x/(4x+1) (s=-1/4), y = x + sigma(x) = 4x^2/(4x+1) = v^2,
   x(v) = (v^2 + v sqrt(1+v^2))/2, and measure the 2-adic slope and the
   archimedean radius in y for:  pure Li_j(-4x),  the inner and outer
   doubly-small companions, and the fold-regular conditional companions.     */
default(parisizemax, 8*10^9);
NTERM = 200;
read("/home/ubuntu/code/math-modular-sources/lattice/catalan_two_classes/00_setup.gp");
GG = Catalan; Z2 = Pi^2/6;
H16 = mkhost(x16,F16); Av = Aof(H16); NMAX = NT-2;
BE = vector(3, j, Bof(H16, mkPhi(Sin ,[[2^(j-1),1]])));
BT = vector(3, j, Bof(H16, mkPhi(Tout,[[2^(j-1),1]])));
rate(sr,n) = abs(polcoeff(sr,n))^(1./n);
v2s(c,k) = if(type(c)!="t_INT" && type(c)!="t_FRAC", "n/a", if(c==0, "zero", valuation(c,2)/(k*1.)));
sig(f) = subst(f, q, -q/(4*q+1) + O(q^(NMAX+1)));
v = 'v;
xofv = (v^2 + v*sqrt(1+v^2+O(v^(2*NMAX))))/2;
toY(g) = subst(g, q, xofv);
inY(sr) = my(r=0); for(k=0,NMAX\2, r += polcoeff(sr,2*k)*q^k); r + O(q^(NMAX\2+1));
{oddchk(sr) = my(bad=-1); for(k=1,NMAX, if(k%2==1 && polcoeff(sr,k)!=0, bad=k;break)); bad;}
{sy(nm, f) =
  my(sp = f + sig(f), SP = toY(sp), YP = inY(SP), K = NMAX\2 - 2);
  print("   ", nm, ":  first odd v-coeff of Sym+ = ", oddchk(SP));
  print("        v2(d_k)/k  k=20,40,",K,": ", [v2s(polcoeff(YP,20),20), v2s(polcoeff(YP,40),40), v2s(polcoeff(YP,K),K)]);
  print("        |d_k|^1/k  k=20,40,",K,": ", [rate(YP,20),rate(YP,40),rate(YP,K)]);
  YP;}
B0in  = BE[1]+4*BE[2]-32*BE[3];
B0out = BT[1]+3*BT[2]-4*BT[3];
Hcond = BE[1] + 7/3*Av;   /* generic rational in place of G: keeps everything exact */                 
Houtc = BT[1]-BT[2] + 11/5*Av;         
Li2 = sum(n=1,NMAX,(-4)^n/n^2*q^n)+O(q^(NMAX+1));
Li1 = sum(n=1,NMAX,(-4)^n/n*q^n)+O(q^(NMAX+1));
print("=== 2-adic slope and radius in y (branch point y=4s=-1, extra point y=-1/2) ===");
sy("A            ", Av);
sy("Li1(-4x)     ", Li1);
sy("Li2(-4x)     ", Li2);
sy("B0in  (1+4V2-32V4)E", B0in);
sy("B0out (1+3V2-4V4)T ", B0out);
print("\n(conditional functions carry an irrational xi; slopes reported for the record)");
sy("Hcond inner  ", Hcond);
sy("Houtc outer  ", Houtc);
quit
