/* Level-16 host: the fold-regular subspace, the target-zero class, 2-adic slopes. */
read("00_setup.gp");
G = Catalan; Z2 = Pi^2/6;
H16 = mkhost(x16,F16); Av = Aof(H16);
NMAX = NT-2;
BE1 = Bof(H16, mkPhi(Sin,[[1,1]]));      /* E     */
BE2 = Bof(H16, mkPhi(Sin,[[2,1]]));      /* V2 E  */
BE4 = Bof(H16, mkPhi(Sin,[[4,1]]));      /* V4 E  */
BT1 = Bof(H16, mkPhi(Tout,[[1,1]]));
BT2 = Bof(H16, mkPhi(Tout,[[2,1]]));
BT4 = Bof(H16, mkPhi(Tout,[[4,1]]));
cmb(c1,c2,c4,U1,U2,U4) = c1*U1+c2*U2+c4*U4;

rate(S,n) = abs(polcoeff(S,n))^(1./n);
{report(nm, Bv, xi) =
  my(Cv = Bv - xi*Av);
  print("\n--- ", nm);
  print("   xi(pred)=",xi);
  print("   b_n/a_n  n=30,60,",NMAX,": ",[polcoeff(Bv,30)/polcoeff(Av,30)*1., polcoeff(Bv,60)/polcoeff(Av,60)*1., polcoeff(Bv,NMAX)/polcoeff(Av,NMAX)*1.]);
  print("   |b_n|^1/n  n=60,100,",NMAX,": ",[rate(Bv,60),rate(Bv,100),rate(Bv,NMAX)]);
  print("   |c_n|^1/n  n=60,100,",NMAX,": ",[rate(Cv,60),rate(Cv,100),rate(Cv,NMAX)]);
  print("   v2(b_n)    n=60,100,",NMAX,": ",[valuation(polcoeff(Bv,60),2),valuation(polcoeff(Bv,100),2),valuation(polcoeff(Bv,NMAX),2)]);
  print("   v2(b_n)/n  n=60,100,",NMAX,": ",[valuation(polcoeff(Bv,60),2)/60.,valuation(polcoeff(Bv,100),2)/100.,valuation(polcoeff(Bv,NMAX),2)/(NMAX*1.)]);
}
print("A(x) = ", Av + O(q^10));
print("|a_n|^1/n n=60,100,",NMAX,": ",[rate(Av,60),rate(Av,100),rate(Av,NMAX)]);
print("v2(a_n) n=60,100: ",[valuation(polcoeff(Av,60),2),valuation(polcoeff(Av,100),2)]);

P2(c1,c2,c4) = c1+c2/4+c4/16;
{ CLS = [ ["E                 ",1,0,0], ["V2 E              ",0,1,0], ["V4 E              ",0,0,1],
          ["Phi_2 = V2E-8V4E  ",0,1,-8], ["Phi_0 = E+4V2E-32V4E",1,4,-32],
          ["E-8V2E (level 8)  ",1,-8,0], ["random 3E+5(V2E-8V4E)",3,5,-40],
          ["generic 1+V2+V4   ",1,1,1] ]; }
for(i=1,#CLS, my(c=CLS[i]); report(c[1], cmb(c[2],c[3],c[4],BE1,BE2,BE4), -P2(c[2],c[3],c[4])*G/2));
print("\n=========== OUTER classes (period P(2)*zeta(2)/2) ===========");
{ CLO = [ ["T                 ",1,0,0], ["V2T               ",0,1,0], ["V4T               ",0,0,1],
          ["T-V2T             ",1,-1,0], ["V2T-8V4T          ",0,1,-8], ["T-5V2T+4V4T       ",1,-5,4] ]; }
for(i=1,#CLO, my(c=CLO[i]); report(c[1], cmb(c[2],c[3],c[4],BT1,BT2,BT4), P2(c[2],c[3],c[4])*Z2/2));

print("\n=========== denominators of Phi_0 companion ===========");
B0 = cmb(1,4,-32,BE1,BE2,BE4);
dn(n)=lcm(vector(n,j,j));
for(k=0,3, my(ok=1,bad=0); for(n=1,NMAX, if(denominator(dn(n)^k*polcoeff(B0,n))!=1, ok=0; if(bad==0,bad=n))); print("   k=",k,": d_n^k B0 integral? ",ok," first failure n=",bad));
print("   B0 first coeffs: ", B0 + O(q^10));
