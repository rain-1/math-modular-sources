/* b7_local.gp -- local exponents of AESZ 207 at 0, the three finite singularities, and infinity. */
default(realprecision, 60);
read("../mum_survey/ops.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
P = OP[4];   /* [P_0(X),...,P_4(X)] */
/* c_k(z) = sum_i p_{i,k} z^i */
ck(k) = sum(i=1,#P, polcoef(P[i],k,'X)*z^(i-1));
/* q_j(z) = z^j sum_k S(k,j) c_k(z) */
S(k,j) = sum(r=0,j, (-1)^(j-r)*binomial(j,r)*r^k)/j!;
{ Q = vector(5, jj, my(j=jj-1); z^j*sum(k=0,4, if(j==0, if(k==0,1,0), S(k,j))*ck(k))); }
print("q_4(z) = ", Q[5]);
print("factor q4/z^4 = ", factor(Q[5]/z^4));
print("q_3(z) = ", Q[4]);
/* indicial at a point z0 (algebraic or rational): work in Q[y]/(minpoly) if needed */
indicial(z0) =
{ my(qs = vector(5,jj, subst(Q[jj], z, z0+t)), mu=10^9, ind=0);
  qs = vector(5,jj, my(s=qs[jj]); if(s==0, 0, s));
  my(ords = vector(5,jj, if(qs[jj]==0, 10^9, valuation(qs[jj], t))));
  for(jj=1,5, if(ords[jj]<10^8, mu = min(mu, ords[jj]-(jj-1))));
  for(jj=1,5, if(ords[jj]-(jj-1)==mu,
     my(lc = polcoef(qs[jj], ords[jj], t), j=jj-1);
     ind += lc*prod(r=0,j-1, rho-r)));
  [mu, ind];
}
print("\n--- z = 1/53248 (double root of R) ---");
{ my(r=indicial(1/53248)); print("  mu=",r[1],"  indicial = ", factor(r[2])); }
print("\n--- z = 0 (MUM) ---");
{ my(r=indicial(0)); print("  mu=",r[1],"  indicial = ", factor(r[2])); }
print("\n--- z = roots of z^2-349/65536 z-1/16777216 (numeric) ---");
{ my(rr=polroots(z^2-349/65536*z-1/16777216));
  for(i=1,2, my(z0=rr[i], q=indicial(z0));
    print("  z0=",z0);
    print("    mu=",q[1],"  indicial roots = ", polroots(q[2]))); }
quit
