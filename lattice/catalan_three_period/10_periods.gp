/* 10_periods.gp -- Eichler period polynomials of the six weight-3 chi_{-4}
   Eisenstein classes on the level-16 host, at every cusp of Gamma_0(16).

   Conventions.  Phi in M_3(Gamma_0(16),chi_{-4}); Theta = D^{-2}(Phi - a_0)
   (positive part only: this is the object the project's companion
   B_Phi = F * Theta uses).  With ahat = Theta + (a_0/2)(2 pi i tau)^2 one has
   D^2 ahat = Phi exactly, hence for g=[a,b;c,d] in Gamma_0(16)

      (c tau+d) ahat(g tau) = chi(d) ( ahat(tau) + alp + bet*tau ).

   Therefore, writing corr(tau) = 2 pi^2 a_0 [ chi(d)^{-1}(c tau+d)(g tau)^2 - tau^2 ],

      R(tau) := chi(d)^{-1} (c tau+d) Theta(g tau) - Theta(tau) = alp + bet*tau + corr(tau).

   corr has a genuine pole at tau=-d/c whenever a_0 != 0 (residue 2 pi^2 a_0/c^2).

   Monodromy of H_xi = F*(Theta - xi) under g:
      H_xi(g tau) - H_xi(tau) = F(tau) * [ R(tau) - xi ( chi(d)(c tau+d) - 1 ) ].
   Fold-regularity at the cusp fixed by g  <=>  R = xi(chi(d)(c tau+d)-1) identically.
*/
default(parisizemax, 8*10^9);
default(realprecision, 130);

NTOP = 6000;
chim(n) = if(n%2==0, 0, if(n%4==1, 1, -1));
ainn = vector(NTOP, n, sumdiv(n, d, chim(n/d)*d^2));   /* inner  E_3^{chi,1}, a_0=0    */
aout = vector(NTOP, n, sumdiv(n, d, chim(d)*d^2));     /* outer  E_3^{1,chi}, a_0=-1/4 */
shft(v, dd) = vector(#v, n, if(n%dd==0, v[n/dd], 0));

CLNAME = ["E      ", "V2E    ", "V4E    ", "T      ", "V2T    ", "V4T    "];
CLVEC  = [ainn, shft(ainn,2), shft(ainn,4), aout, shft(aout,2), shft(aout,4)];
CLA0   = [0, 0, 0, -1/4, -1/4, -1/4];

/* Theta and its tau-derivative at a point of the upper half plane */
thv(v, tau, NN) = my(qq=exp(2*Pi*I*tau), s=0., p=1.); for(n=1,NN, p*=qq; if(v[n]!=0, s += v[n]/n^2*p)); s;
thd(v, tau, NN) = my(qq=exp(2*Pi*I*tau), s=0., p=1.); for(n=1,NN, p*=qq; if(v[n]!=0, s += v[n]/n*p));  2*Pi*I*s;
PRECDIG = 130;
nterms(tau) = my(y=imag(tau)); ceil((PRECDIG+25)*log(10)/(2*Pi*y));

/* the six cusps of Gamma_0(16) with a generator of the stabiliser */
/* [name, width, gamma] */
{CUSPS = [
  ["0  ", 16, [1,0;-16,1]],
  ["1/2", 4,  [-7,4;-16,9]],
  ["1/4", 1,  [-3,1;-16,5]],
  ["3/4", 1,  [-11,9;-16,13]],
  ["1/8", 1,  [-7,1;-64,9]]
];}

/* corr(tau) and corr'(tau) */
{corrf(gm, a0, tau) = my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2], ch=chim(d));
  2*Pi^2*a0*( (a*tau+b)^2/(ch*(c*tau+d)) - tau^2 );}
{corrd(gm, a0, tau) = my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2], ch=chim(d));
  2*Pi^2*a0*( (2*a*(a*tau+b)*(c*tau+d) - c*(a*tau+b)^2)/(ch*(c*tau+d)^2) - 2*tau );}

/* period polynomial (alp,bet) of class j at the cusp fixed by g */
{perpol(gm, v, a0, frac) =
  my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2], ch=chim(d));
  my(p0 = (a-d)/(2*c));                 /* fixed point */
  my(w0 = -c/2 - I*abs(c)*frac);        /* Im(w)<0  <=>  Im(tau)>0 */
  my(t0 = p0 + 1/w0, t1 = (a*t0+b)/(c*t0+d));
  my(N0 = nterms(t0), N1 = nterms(t1));
  my(NN = max(N0,N1));
  if(NN > NTOP, error("need more terms: ", NN));
  my(R  = (c*t0+d)*thv(v,t1,NN)/ch - thv(v,t0,NN));
  my(Rp = ( c*thv(v,t1,NN) + thd(v,t1,NN)/(c*t0+d) )/ch - thd(v,t0,NN));
  my(bet = Rp - corrd(gm,a0,t0));
  my(alp = R - corrf(gm,a0,t0) - bet*t0);
  [alp, bet, t0];
}

print("=== period polynomials (alp,bet) at each cusp of Gamma_0(16) ===");
print("Theta = D^{-2}(Phi-a_0);  R(tau)-corr(tau) = alp + bet*tau");
GG = Catalan; Z2 = Pi^2/6;
{for(k=1,#CUSPS,
  my(cu=CUSPS[k], gm=cu[3], c=gm[2,1], d=gm[2,2]);
  print("\n--- cusp ", cu[1], "  width ", cu[2], "  gamma=", gm, "  chi(d)=", chim(d));
  print("    fold-regularity condition:  ", c, "*alp - ", d-1, "*bet = 0 ;   xi = bet/", c);
  for(j=1,6,
    my(r = perpol(gm, CLVEC[j], CLA0[j], 0.5));
    my(r2 = perpol(gm, CLVEC[j], CLA0[j], 0.75));
    print("  ", CLNAME[j], " alp=", r[1], "  bet=", r[2]);
    print("           check (2nd pt): dalp=", abs(r[1]-r2[1]), " dbet=", abs(r[2]-r2[2]));
  );
);}
quit
