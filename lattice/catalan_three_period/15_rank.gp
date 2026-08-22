/* 15_rank.gp -- Q(y)-linear independence of the small functions and of their
   CDT-shaped orbits on the level-16 host.

   CDT's 7 conditional functions (verified in lattice/cdt_finder/indep_check2.py)
   are the orbit of ONE symmetrised conditional function g = Sym^+ H :
        g, g', g'', g''' ,  int g dy,  int g dy/y,  int g dy/y^2 .
   Here we build the same 7-member orbit for each of the three small functions
   of the level-16 host (inner doubly-small, outer doubly-small, conditional)
   and test Q(y)-independence by the rank of {y^j f_i} mod a large prime.     */
default(parisizemax, 8*10^9);
NTERM = 200;
read("/home/ubuntu/code/math-modular-sources/lattice/catalan_two_classes/00_setup.gp");
H16 = mkhost(x16,F16); Av = Aof(H16); NMAX = NT-2;
BE = vector(3, j, Bof(H16, mkPhi(Sin ,[[2^(j-1),1]])));
BT = vector(3, j, Bof(H16, mkPhi(Tout,[[2^(j-1),1]])));
sig(f) = subst(f, q, -q/(4*q+1) + O(q^(NMAX+1)));
v = 'v;
xofv = (v^2 + v*sqrt(1+v^2+O(v^(2*NMAX))))/2;
KY = NMAX\2 - 2;
{symp(f) = my(sp=f+sig(f), SP=subst(sp,q,xofv), r=0);
  for(k=0,KY, r += polcoeff(SP,2*k)*q^k); r + O(q^(KY+1));}
{symm(f) = my(sm=f-sig(f), od=q-sig(q));            /* odd part / (x - sigma x) */
  my(sr = sm/od, SR=subst(sr,q,xofv), r=0);
  for(k=0,KY, r += polcoeff(SR,2*k)*q^k); r + O(q^(KY+1));}
KZ = KY - 8;
trun(f) = my(r=0); for(k=0,KZ, r += polcoeff(f,k)*q^k); r + O(q^(KZ+1));
dy(f) = my(r=0); for(k=1,KZ, r += k*polcoeff(f,k)*q^(k-1)); r + O(q^(KZ+1));
iy(f) = my(r=0); for(k=0,KZ-1, r += polcoeff(f,k)/(k+1)*q^(k+1)); r + O(q^(KZ+1));
{iyw(f,w) = my(r=0);      /* int f dy / y^w  (needs coefficients 0..w-1 to vanish) */
  for(k=w,KZ-1, r += polcoeff(f,k)/(k+1-w)*q^(k+1-w)); r + O(q^(KZ+1));}
GAM = 7/3;                                   /* generic rational standing for G   */
ZET = 11/5;                                  /* generic rational standing for z(2) */
B0in  = BE[1]+4*BE[2]-32*BE[3];
B0out = BT[1]+3*BT[2]-4*BT[3];
Hin   = BE[1] + GAM/2*Av;                    /* inner conditional  (xi=-G/2)      */
Li2 = sum(n=1,NMAX,(-4)^n/n^2*q^n)+O(q^(NMAX+1));
Li3 = sum(n=1,NMAX,(-4)^n/n^3*q^n)+O(q^(NMAX+1));
Li1 = sum(n=1,NMAX,(-4)^n/n*q^n)+O(q^(NMAX+1));
{orb(f) = my(g=trun(f));   /* CDT's 7-member orbit */
  [g, dy(g), dy(dy(g)), dy(dy(dy(g))), iy(g), iyw(g,1), iyw(g,2)];}
g1 = symp(B0in); g2 = symp(B0out); g3 = symp(Hin);
print("Sym+ series check: g1,g2,g3 first coeffs");
print("  g1 = ", g1 + O(q^6)); print("  g2 = ", g2 + O(q^6)); print("  g3 = ", g3 + O(q^6));
O1 = orb(g1); O2 = orb(g2); O3 = orb(g3);
PU = [trun(1+O(q^(KY+1))), trun(symp(Li2)), trun(symm(Li2)), trun(symp(Li3)), trun(symm(Li3)), trun(symp(Li1)), trun(symm(Li1))];
PP = 2^61-1;
{rk(fl, D) = my(rows=[], nc=KZ-1);
  for(i=1,#fl, for(j=0,D, my(f=fl[i], r=vector(nc));
    for(k=0,nc-1, my(kk=k-j); r[k+1] = if(kk>=0, Mod(polcoeff(f,kk),PP), Mod(0,PP)));
    rows = concat(rows,[r])));
  matrank(Mat(matconcat(rows~)~));}
print("\n=== rank of {y^j f_i} mod 2^61-1 (need rank = #f * (D+1)) ===");
{for(D=0,3,
  print("  D=",D,":  small 21: ", rk(concat(concat(O1,O2),O3), D), " / ", 21*(D+1),
        "   |  small+pure 28: ", rk(concat(concat(concat(O1,O2),O3),PU), D), " / ", 28*(D+1),
        "   |  pure 7: ", rk(PU,D), " / ", 7*(D+1));)}
print("\n=== rank of the three generators alone (and with A) ===");
{for(D=0,5, print("  D=",D,": {g1,g2,g3}: ", rk([trun(g1),trun(g2),trun(g3)],D), "/", 3*(D+1),
                  "   {g1,g2,g3,Sym+A}: ", rk([trun(g1),trun(g2),trun(g3),trun(symp(Av))],D), "/", 4*(D+1)));}
quit
