/* broad (order,degree) scan for the level-12 hosts: what order does A(t) satisfy? */
default(parisizemax, 24000000000);
default(threadsizemax, 8000000000);
\p 40
if(type(MM)!="t_INT", MM = 420);
M = MM;
q = 'q; ee(d) = eta(q^d + O(q^(M+2)));
PP = 2305843009213693951;
mkrow(cv, dv, T) =
{ my(g, PHI, DT, F, QT, A, s);
  g = vector(M);
  for(n=1,M, s=0; for(i=1,#dv, if(n%dv[i]==0, s += cv[i]*sigma(n/dv[i],5))); g[n]=s);
  PHI = sum(n=1,M, g[n]*q^n) + O(q^(M+1));
  DT = q*deriv(T,q); F = PHI/DT; QT = serreverse(T);
  A = subst(F,q,QT);
  vector(M, i, polcoeff(A, i-1));
}
pfrow(V, R, D, n) =
{ my(L=List(), m, co);
  for(k=0,R, for(j=0,D,
    m = n - j + k;
    if(m < 0 || m > #V-1, co = 0, co = prod(i=0,k-1, n-j+k-i) * V[m+1]);
    listput(L, co)));
  Vec(L);
}
bigscan(name, V) =
{ my(NN=#V, Vm, found=0);
  Vm = vector(NN, i, Mod(V[i], PP));
  print("=== ", name, "  (terms to n=", NN-1, ")");
  for(R=2, 12,
    for(D=1, 90,
      my(nv=(R+1)*(D+1));
      if(nv+30 > NN, break());
      my(nhi=min(NN-1, nv+60),
         K = #matker(matconcat(vector(nhi+1, i, pfrow(Vm,R,D,i-1))~)));
      if(K>0, print("   FIT: order ", R, " degree ", D, "  kernel dim ", K,
                    "  [eqns ", nhi+1, ", unknowns ", nv, "]"); found=1; break())));
  if(!found, print("   no fit with order<=12 in the searched range"));
}
WW  = q*(ee(1)*ee(12)/(ee(3)*ee(4)))^4 + O(q^(M+1));  TD = WW/(1+WW)^2;
HH  = (ee(1)^3*ee(4)*ee(6)^2/(ee(2)^2*ee(3)*ee(12)^3))/q + O(q^(M+1)); TX = HH/((HH+3)*(HH+4));
XX  = q*ee(2)*ee(16)^2/(ee(1)^2*ee(8)) + O(q^(M+1)); T16 = XX/(8*XX^2+2*XX+1);
c12 = [1,-113,567,112,-1863,1296]; d12 = [1,2,3,4,6,12];
c16 = [1,-85,1428,-5440,4096];     d16 = [1,2,4,8,16];
bigscan("level-12 Domb", mkrow(c12,d12,TD));
bigscan("level-12 h_12", mkrow(c12,d12,TX));
bigscan("level-16     ", mkrow(c16,d16,T16));
print("DONE"); quit;
