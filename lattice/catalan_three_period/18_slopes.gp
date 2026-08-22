/* 18_slopes.gp -- the 2-adic slope in y of EVERY object entering the inventory,
   including the host A and the raw companions.  The point: y = 4x^2/(4x+1) is
   NOT an integral coordinate (x(v) = (v^2+v sqrt(1+v^2))/2 has 2-power
   denominators), so Sym^+ of an integral series need not be integral in y.     */
default(parisizemax, 8*10^9);
NTERM = 200;
read("/home/ubuntu/code/math-modular-sources/lattice/catalan_two_classes/00_setup.gp");
H16 = mkhost(x16,F16); Av = Aof(H16); NMAX = NT-2;
BE = vector(3, j, Bof(H16, mkPhi(Sin ,[[2^(j-1),1]])));
BT = vector(3, j, Bof(H16, mkPhi(Tout,[[2^(j-1),1]])));
sig(f) = subst(f, q, -q/(4*q+1) + O(q^(NMAX+1)));
v = 'v; xofv = (v^2 + v*sqrt(1+v^2+O(v^(2*NMAX))))/2;
{symp(f) = my(sp=f+sig(f), SP=subst(sp,q,xofv), K=NMAX\2-2, r=0);
  for(k=0,K, r += polcoeff(SP,2*k)*q^k); r + O(q^(K+1));}
rate(sr,n) = abs(polcoeff(sr,n))^(1./n);
v2s(c,k) = if(c==0, "zero", valuation(c,2)/(k*1.));
K = NMAX\2-3;
{rp(nm, f) = my(Y=symp(f));
  print("  ", nm, "  v2(d_k)/k at k=20,40,",K,": ", [v2s(polcoeff(Y,20),20), v2s(polcoeff(Y,40),40), v2s(polcoeff(Y,K),K)],
        "   |d_k|^1/k: ", [rate(Y,40), rate(Y,K)]);}
Li1 = sum(n=1,NMAX,(-4)^n/n*q^n)+O(q^(NMAX+1));
Li2 = sum(n=1,NMAX,(-4)^n/n^2*q^n)+O(q^(NMAX+1));
Li3 = sum(n=1,NMAX,(-4)^n/n^3*q^n)+O(q^(NMAX+1));
print("=== 2-adic slopes in y = 4x^2/(4x+1) ===");
rp("A        (host)          ", Av);
rp("1                        ", 1+O(q^(NMAX+1)));
rp("x                        ", q+O(q^(NMAX+1)));
rp("Li2(-4x)                 ", Li2);
rp("Li3(-4x)                 ", Li3);
rp("B_E   (raw inner cmpn)   ", BE[1]);
rp("B_T-B_V2T (raw outer)    ", BT[1]-BT[2]);
rp("H_cond = B_E + (7/3)A    ", BE[1]+7/3*Av);
rp("H_cond = B_E + (1/8)A    ", BE[1]+1/8*Av);
rp("B0in  = (1+4V2-32V4)E    ", BE[1]+4*BE[2]-32*BE[3]);
rp("B0out = (1+3V2-4V4)T     ", BT[1]+3*BT[2]-4*BT[3]);
quit
