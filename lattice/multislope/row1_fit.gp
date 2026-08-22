/* Row 1: A_n = [t^n] F(t)^3, F = sum e_n t^n (Zagier E, level 6-ish).
   (a) minimal recurrence scan + exact reconstruction
   (b) characteristic roots                                            */
default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_fit.log";
W(s) = write(LOG, s);

NMAX = 420;
/* Zagier E */
ev = vector(NMAX+2);
ev[1] = 1; ev[2] = 4;
for(n=1, NMAX, ev[n+2] = ((12*n^2+12*n+4)*ev[n+1] - 32*n^2*ev[n])/(n+1)^2);
W("=== Row 1: Zagier E and A_n = [t^n] F^3 ===");
W(Str("e_0..e_10 = ", vector(11,i,ev[i])));
bad = 0; for(i=1,NMAX+2, if(denominator(ev[i])!=1, bad++));
W(Str("all e_n integral over n<=", NMAX+1, "? ", bad==0));

/* A = coefficients of F^3 */
FS = sum(n=0, NMAX, ev[n+1]*'t^n) + O('t^(NMAX+1));
AS = FS^3;
Av = vector(NMAX+1, i, polcoeff(AS, i-1));
W(Str("A_0..A_10 = ", vector(11,i,Av[i])));

/* ---- mod-p kernel-dimension scan: sum_{i=0}^{r} sum_{j=0}^{DG} c_{i,j} n^j A_{n-i} = 0 ---- */
pp = 2^61 - 1;
Ap = vector(NMAX+1, i, Mod(Av[i], pp));
W("");
W("--- kernel dimensions mod p=2^61-1 for (order r, poly degree DG) ---");
W("   r \\ DG :  2   3   4   5   6   7   8");
for(r=2, 8, \
  row = Str("   r=", r, "   : "); \
  for(dg=2, 8, \
    nun = (r+1)*(dg+1); \
    nrow = nun + 40; \
    if(r + nrow <= NMAX, \
      mt = matrix(nrow, nun, a, b, \
            my(i = (b-1)\(dg+1), j = (b-1)%(dg+1), n = r + a); \
            Mod(n,pp)^j * Ap[n-i+1]); \
      kd = #matker(mt), kd = -1); \
    row = Str(row, kd, "   ")); \
  W(row));
W("");
W("DONE-SCAN");
quit;
