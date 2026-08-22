/* 13_outer_dbl.gp -- the OUTER fold-regular plane on the level-16 host and its
   target-zero (doubly small) direction.

   From 12_exact.gp (period polynomials at the fold cusp 1/2):
     outer classes: alp = 4 zeta(2) Q,  bet = -8 zeta(2) Q,  Q = P(0)+3P(1)-3P(2)
                    a_0(Phi) = -P(0)/4
     => the companion B_Phi is fold-regular  <=>  P(0)=0,
        and then  xi_inf = -bet/16 = (zeta(2)/2)(P(0)+3P(1)-3P(2)) = (3/2)zeta(2)(P(1)-P(2)).
     Target zero  <=>  P(0)=0 and P(1)=P(2)  =>  Phi = (1+3V2-4V4)T  (up to scale).
   This script tests those predictions numerically on the actual power series.  */
default(parisizemax, 8*10^9);
NTERM = 220;
read("/home/ubuntu/code/math-modular-sources/lattice/catalan_two_classes/00_setup.gp");
GG = Catalan; Z2 = Pi^2/6;
H16 = mkhost(x16,F16); Av = Aof(H16); NMAX = NT-2;
BT = vector(3, j, Bof(H16, mkPhi(Tout,[[2^(j-1),1]])));   /* T, V2T, V4T */
BE = vector(3, j, Bof(H16, mkPhi(Sin ,[[2^(j-1),1]])));   /* E, V2E, V4E */
rate(sr,n) = abs(polcoeff(sr,n))^(1./n);
Pv(c,s) = c[1] + c[2]*2^(-s) + c[3]*4^(-s);
xiout(c) = Z2/2*(Pv(c,0)+3*Pv(c,1)-3*Pv(c,2));
xiinn(c) = -GG/2*Pv(c,2);
dn(n) = lcm(vector(n,j,j));
{rep(nm, Bv, xi) = my(Cv = Bv - xi*Av);
  print("\n--- ", nm, "   xi = ", xi);
  print("   b_n/a_n   n=60,120,",NMAX,": ",[polcoeff(Bv,60)/polcoeff(Av,60)*1., polcoeff(Bv,120)/polcoeff(Av,120)*1., polcoeff(Bv,NMAX)/polcoeff(Av,NMAX)*1.]);
  print("   |b_n|^1/n n=60,120,",NMAX,": ",[rate(Bv,60),rate(Bv,120),rate(Bv,NMAX)]);
  print("   |c_n|^1/n n=60,120,",NMAX,": ",[rate(Cv,60),rate(Cv,120),rate(Cv,NMAX)],"   (fold 4, next 2sqrt2=",2*sqrt(2),")");
  print("   v2(b_n)   n=60,120,",NMAX,": ",[valuation(polcoeff(Bv,60),2),valuation(polcoeff(Bv,120),2),valuation(polcoeff(Bv,NMAX),2)]);
  print("   v2(b_n)/n            : ",[valuation(polcoeff(Bv,60),2)/60.,valuation(polcoeff(Bv,120),2)/120.,valuation(polcoeff(Bv,NMAX),2)/(NMAX*1.)]);
}
print("host: |a_n|^1/n n=60,120,",NMAX,": ",[rate(Av,60),rate(Av,120),rate(Av,NMAX)]);
print("\n############ OUTER classes, fold-regular xi = (3/2)zeta(2)(P(1)-P(2)) ############");
{CLO = [ ["T            ",[1,0,0]], ["V2T          ",[0,1,0]], ["V4T          ",[0,0,1]],
         ["T-V2T        ",[1,-1,0]], ["T-5V2T+4V4T  ",[1,-5,4]], ["V2T-8V4T     ",[0,1,-8]],
         ["OUT0=T+3V2T-4V4T",[1,3,-4]], ["2T+V2T-3V4T  ",[2,1,-3]], ["T+V2T-2V4T   ",[1,1,-2]] ];}
{for(i=1,#CLO, my(c=CLO[i][2], Bv=sum(j=1,3,c[j]*BT[j]));
  print("\n   P(0)=",Pv(c,0)," P(1)=",Pv(c,1)," P(2)=",Pv(c,2)," a0=",-Pv(c,0)/4);
  rep(CLO[i][1], Bv, xiout(c)));}
print("\n############ the two target-zero classes ############");
B0in  = BE[1]+4*BE[2]-32*BE[3];
B0out = BT[1]+3*BT[2]-4*BT[3];
rep("INNER Phi0 = (1+4V2-32V4)E ", B0in, 0);
rep("OUTER Phi0 = (1+3V2-4V4)T  ", B0out, 0);
print("\n   B0out first coeffs: ", B0out + O(q^12));
{for(k=0,3, my(ok=1,bad=0); for(n=1,NMAX, if(denominator(dn(n)^k*polcoeff(B0out,n))!=1, ok=0; if(bad==0,bad=n)));
   print("   d_n^",k," B0out integral? ",ok,"  first failure n=",bad));}
{for(k=0,3, my(ok=1,bad=0); for(n=1,NMAX, if(denominator(dn(n)^k*polcoeff(B0in,n))!=1, ok=0; if(bad==0,bad=n)));
   print("   d_n^",k," B0in  integral? ",ok,"  first failure n=",bad));}
quit
