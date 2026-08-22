/* Common setup: hosts, weight-3 chi_{-4} Eisenstein classes, companions. */
q = 'q;
if(type(NTERM)!="t_INT", NTERM = 120);
NT = NTERM;

et(k) = eta(q^k + O(q^(NT+2)));

/* ---- hosts (eta quotients) ---- */
t8  = q*et(1)^4*et(4)^2*et(8)^4/et(2)^10  + O(q^NT);
F8  = et(2)^10/(et(1)^4*et(4)^4)          + O(q^NT);
x16 = q*et(2)*et(16)^2/(et(1)^2*et(8))    + O(q^NT);
F16 = et(2)^4/et(4)^2                     + O(q^NT);

/* ---- weight-3 chi_{-4} Eisenstein series (coefficient vectors, n>=1) ---- */
chi(n) = if(n%2==0,0, if(n%4==1,1,-1));
Sin  = vector(NT-1, m, sumdiv(m,d, chi(m/d)*d^2));   /* inner: psi=chi_{-4}, phi=1  */
Tout = vector(NT-1, m, sumdiv(m,d, chi(d)*d^2));     /* outer: psi=1, phi=chi_{-4}  */
T0const = -1/2;   /* L(-2,chi_{-4}) = -1/2 : constant term of the outer series */

Vd(v,d) = vector(#v, m, if(m%d==0, v[m/d], 0));
mkPhi(base, pl) = sum(i=1,#pl, pl[i][2]*Vd(base,pl[i][1]));
Peval(pl,s) = sum(i=1,#pl, pl[i][2]*pl[i][1]^(-s));

ser(v)   = sum(m=1,#v, v[m]*q^m) + O(q^(#v+1));
Theta(v) = sum(m=1,#v, v[m]/m^2*q^m) + O(q^(#v+1));   /* D^{-2} on positive part */

/* host = [t(q), F(q), q(t)] */
mkhost(tq,Fq) = [tq, Fq, serreverse(tq)];
Aof(H)    = subst(H[2], q, H[3]);
Bof(H,v)  = subst(H[2]*Theta(v), q, H[3]);
