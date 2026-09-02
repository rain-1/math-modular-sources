/* lib.gp -- cusp periods of the Eichler integrals of the Apery-row Eisenstein
   sources.  Conventions:
     E_k^{psi,phi} = c_0 + sum_m (sum_{e|m} psi(m/e) phi(e) e^{k-1}) q^m ,  k = r+1
     Phi = sum_d c_d E_{r+1}^{psi,phi}(d tau)                (oldform combination)
     Theta = D^{-r} Phi = sum_{m>=1} c(m) m^{-r} q^m ,       D = q d/dq
   At the cusp a/c put q = zeta e^{-t}, zeta = e(a/c), t -> 0+.  Then
     Theta = rho/t + Pi_{a/c} + O(t log t).
   Grouping D^{-r}E = sum_f psi(f) f^{-r} Lam_phi(q^f),  Lam_phi(z) = N_phi(z)/(1-z^Q),
   N_phi(z) = sum_{a=1}^{Q} phi(a) z^a, Q = period of phi, gives
     Pi   = sum_d c_d d^{-r}   sum_f psi(f) f^{-r}   w_phi(zeta^{df})
     rho  = sum_d c_d d^{-r-1} sum_f psi(f) f^{-r-1} n_phi(zeta^{df})
   with, at z0 with z0^Q = 1:  w = N(z0)/2 - z0 N'(z0)/Q ,  n = N(z0)/Q ;
        otherwise             w = N(z0)/(1-z0^Q) ,          n = 0 .
   The f-sums are periodic mod M = lcm(c/gcd(c,d), period(psi)) and are evaluated
   exactly by Hurwitz zeta:  sum_f g(f) f^{-s} = M^{-s} sum_{j=1}^M g(j) zetahurwitz(s,j/M).
                                                                             */

/* ---- characters / periodic functions: ch = [Q, [v(1),...,v(Q)]] --------- */
pv(ch, j) = ch[2][((j-1)%ch[1])+1];
Nph(ch, z) = sum(a=1, ch[1], ch[2][a]*z^a);
Ndp(ch, z) = sum(a=1, ch[1], a*ch[2][a]*z^(a-1));
EPS = 10^(-(default(realprecision)\2));
wph(ch, z0) = if(abs(z0^ch[1]-1) < EPS, Nph(ch,z0)/2 - z0*Ndp(ch,z0)/ch[1], Nph(ch,z0)/(1-z0^ch[1]));
nph(ch, z0) = if(abs(z0^ch[1]-1) < EPS, Nph(ch,z0)/ch[1], 0);

ee(x) = exp(2*Pi*I*x);

/* ---- one Eisenstein block E^{psi,phi}(d tau) at cusp a/c ---------------- */
/* returns [Pi_block, rho_block] for the block with d = 1 and twist alpha = a/c */
blk(psi, phi, r, a, c) = {
  my(g = gcd(a,c), aa = a/g, cc = c/g, M = lcm(cc, psi[1]), SP = 0, SR = 0, z);
  for(j = 1, M,
    z = ee((aa*j)/cc);
    SP += pv(psi,j) * wph(phi, z) * zetahurwitz(r,   j/M);
    SR += pv(psi,j) * nph(phi, z) * zetahurwitz(r+1, j/M));
  [SP/M^r, SR/M^(r+1)];
}

/* ---- a source: src = [name, r, N, psi, phi, dlist, clist] --------------- */
cperiod(src, a, c) = {
  my(r=src[2], psi=src[4], phi=src[5], dl=src[6], cl=src[7], P=0, R=0, b);
  for(i=1, #dl, b = blk(psi, phi, r, a*dl[i], c); P += cl[i]*dl[i]^(-r)*b[1]; R += cl[i]*dl[i]^(-r-1)*b[2]);
  [P, R];
}

/* ---- cusp representatives of Gamma_0(n): list of [a,c], c|n ------------- */
cusplist(n) = {my(res=List()); fordiv(n, c, my(g=gcd(c,n/c), seen=List()); for(a=1, max(g,1), if(gcd(a,c)==1, my(new=1); for(k=1,#seen, if((a-seen[k])%g==0, new=0; break)); if(new, listput(seen,a); listput(res,[a,c]))))); Vec(res);}

/* ---- the twelve table sources + the Gamma_1(5) newcomers ---------------- */
TRIV  = [1, [1]];
CHM3  = [3, [1,-1,0]];
CHM4  = [4, [1,0,-1,0]];
CH5   = [5, [1,-1,-1,1,0]];
WD    = [5, [1,-2,2,-1,0]];        /* re(psi_4) - 2 im(psi_4) */
R1v   = [5, [2,0,0,-2,0]];         /* 2 re(psi_4)   */
R2v   = [5, [0,-2,2,0,0]];         /* -2 im(psi_4)  */
PS4   = [5, [1,I,-I,-1,0]];        /* quartic psi_4, psi_4(2) = i */
PS4B  = [5, [1,-I,I,-1,0]];        /* its conjugate */

{SRC = [
 ["A",     2,  6, TRIV, CHM3, [1,2],       [1,-1]],
 ["B",     2, 36, CHM3, TRIV, [1,2,4],     [1,-6,-8]],
 ["C",     2,  6, CHM3, TRIV, [1,2],       [1,-8]],
 ["D",     2,  5, TRIV, WD,   [1],         [1]],
 ["E",     2,  8, CHM4, TRIV, [1,2],       [1,-8]],
 ["F",     2, 12, CHM3, TRIV, [1,2,4],     [1,-7,-8]],
 ["alpha", 3, 12, TRIV, TRIV, [1,2,3,4,6,12], [1,-17,-9,16,153,-144]],
 ["gamma", 3,  6, TRIV, TRIV, [1,2,3,6],   [1,-28,63,-36]],
 ["delta", 3, 12, TRIV, TRIV, [1,2,3,4,6,12], [1,-14,-1,16,14,-16]],
 ["eps",   3,  8, TRIV, TRIV, [1,2,4,8],   [1,-21,84,-64]],
 ["zeta",  3,  9, CHM3, CHM3, [1],         [1]],
 ["eta",   3, 20, CH5,  TRIV, [1,2,4],     [1,-14,-16]]
];}
/* extra Gamma_1(5) weight-three directions (inner quartic) */
{XTRA = [
 ["R1",    2,  5, TRIV, R1v,  [1], [1]],
 ["R2",    2,  5, TRIV, R2v,  [1], [1]],
 ["E3ps",  2,  5, PS4,  TRIV, [1], [1]],
 ["E3psb", 2,  5, PS4B, TRIV, [1], [1]]
];}

srcbyname(nm) = {my(L=concat(SRC,XTRA)); for(i=1,#L, if(L[i][1]==nm, return(L[i]))); error("no source ",nm);}

/* ---- direct Abel check: Theta(zeta e^{-t}) by summation ---------------- */
cmc(src, m) = {my(r=src[2], psi=src[4], phi=src[5], dl=src[6], cl=src[7], s=0);
  for(i=1,#dl, if(m % dl[i] == 0, my(mm = m/dl[i]);
    s += cl[i]*sumdiv(mm, e, pv(psi, mm/e)*pv(phi, e)*e^r))); s;}
abel(src, a, c, tt, MX) = {my(r=src[2], z=ee(a/c), s=0); for(m=1, MX, s += cmc(src,m)/m^r*z^m*exp(-m*tt)); s;}
