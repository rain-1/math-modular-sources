\\ eta-quotient q-expansions and Apery-system utilities
PREC = if(type(PRECSET)=="t_INT",PRECSET,220);
ser(v) = v + O(q^PREC);
\\ E(d) = eta(d tau)/q^(d/24) as a power series in q
E(d) = eta(q^d + O(q^PREC));
\\ etaq(list of [d,e]) : full eta quotient including q-power prefix
etaq(L) = { my(s=0, r=1+O(q^PREC));
  for(i=1,#L, s += L[i][1]*L[i][2]; r *= E(L[i][1])^L[i][2]);
  if(s%24, error("noninteger q power ",s/24));
  q^(s/24)*r };
Dq(f) = q*f';
\\ inverse theta on zero-constant-term series
Dqinv(f, k) = { my(r=O(q^PREC)); for(m=1,PREC-1, r += polcoeff(f,m)/m^k*q^m); r };
Vd(f,d) = subst(f, q, q^d) + O(q^PREC);
