default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_fit2.log";
W(s) = write(LOG, s);
NMAX = 420;
ev = vector(NMAX+2); ev[1]=1; ev[2]=4;
for(n=1, NMAX, ev[n+2] = ((12*n^2+12*n+4)*ev[n+1] - 32*n^2*ev[n])/(n+1)^2);
FS = sum(n=0, NMAX, ev[n+1]*'t^n) + O('t^(NMAX+1));
AS = FS^3;
Av = vector(NMAX+1, i, polcoeff(AS, i-1));
pp = 2^61 - 1; Ap = vector(NMAX+1, i, Mod(Av[i], pp));
W("--- extended kernel-dim scan mod p, DG up to 16 ---");
for(r=2, 7, row = Str("   r=", r, " : "); \
  for(dg=2, 16, nun=(r+1)*(dg+1); nrow=nun+40; \
    if(r+nrow<=NMAX, mt = matrix(nrow, nun, a, b, my(i=(b-1)\(dg+1), j=(b-1)%(dg+1), n=r+a); Mod(n,pp)^j*Ap[n-i+1]); kd=#matker(mt), kd=-1); \
    row=Str(row, kd, " ")); W(row));

/* ---- exact reconstruction at r=6, DG=4 ---- */
r = 6; dg = 4; nun = (r+1)*(dg+1); nrow = nun + 25;
mt = matrix(nrow, nun, a, b, my(i=(b-1)\(dg+1), j=(b-1)%(dg+1), n=r+a); n^j*Av[n-i+1]);
kk = matker(mt);
W(""); W(Str("exact kernel dim at (r,DG)=(6,4): ", #kk));
vv = kk[,1];
den = 1; for(i=1,#vv, den = lcm(den, denominator(vv[i])));
vv = vv*den; g = 0; for(i=1,#vv, g = gcd(g, vv[i])); if(g!=0, vv = vv/g);
if(vv[1] < 0, vv = -vv);
Qp = vector(r+1, i, sum(j=0, dg, vv[(i-1)*(dg+1)+j+1]*'n^j));
W("--- the minimal recurrence:  sum_{i=0}^{6} Q_i(n) A_{n-i} = 0 ---");
for(i=0,r, W(Str("  Q_", i, "(n) = ", Qp[i+1], "     factored: ", factor(Qp[i+1]))));
/* verify on all available terms */
okv = 1;
for(n=r, NMAX, s = sum(i=0,r, subst(Qp[i+1],'n,n)*Av[n-i+1]); if(s!=0, okv=0; W(Str("  FAIL at n=",n))));
W(Str("verified for n=6..", NMAX, ": ", okv));
/* characteristic polynomial: leading coeff in n */
cp = sum(i=0, r, polcoeff(Qp[i+1], dg, 'n) * 'x^(r-i));
W(""); W(Str("characteristic polynomial (n->oo): ", cp));
W(Str("  factored: ", factor(cp)));
W(Str("  roots: ", polroots(cp)));
W(Str("  trailing/leading: ", polcoeff(cp,0,'x)/polcoeff(cp,r,'x)));
/* also low-order coefficients for kappa etc */
W(""); W(Str("Q_0 = ", factor(Qp[1])));
W(Str("Q_6 = ", factor(Qp[7])));
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_rec.txt", Str("QROW1 = ", Qp));
W("DONE");
quit;
