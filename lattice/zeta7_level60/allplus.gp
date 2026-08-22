\\ Which purified weight-8 level-60 sources live on the ONLY genus-zero quotient,
\\ X_0(60)* = X_0(60)/<W_3,W_4,W_5> ?  Need W_Q Phi = +Phi for all Q | 60.
dl = divisors(60);
P(c,s) = sum(i=1,#dl, c[i]*dl[i]^(-s));
pure(c) = vector(5,j,P(c,2*j-2))==[0,0,0,0,0];
{AL(c,Q) = my(k=8, out=vector(#dl));
  for(i=1,#dl, my(d=dl[i], a=gcd(d,Q^10), dn, f, j=0);
    dn = d*Q/a^2; f = (Q/a^2)^(k/2);
    for(t=1,#dl, if(dl[t]==dn, j=t)); out[j] = out[j] + f*c[i]);
  out;}
{ALM(Q) = matrix(#dl,#dl, r,c, 0);}
{ALmat(Q) = my(A=matrix(#dl,#dl));
  for(i=1,#dl, my(e=vector(#dl)); e[i]=1; my(w=AL(e,Q));
    for(j=1,#dl, A[j,i]=w[j])); A;}
M = matrix(5,#dl, r,i, dl[i]^(-(2*(r-1))));
A3=ALmat(3); A4=ALmat(4); A5=ALmat(5); Id=matid(#dl);
print("dim(oldform, all W_Q = +1)            = ", #matker(matconcat([A3-Id;A4-Id;A5-Id])));
print("dim(purified, all W_Q = +1)           = ", #matker(matconcat([M;A3-Id;A4-Id;A5-Id])));
{for(sgn3=-1,1, for(sgn4=-1,1, for(sgn5=-1,1,
  if(sgn3!=0 && sgn4!=0 && sgn5!=0,
   my(K=matker(matconcat([M;A3-sgn3*Id;A4-sgn4*Id;A5-sgn5*Id])));
   print("  signs (W3,W4,W5)=(",sgn3,",",sgn4,",",sgn5,") -> dim purified = ",#K,
     if(#K>0, concat(["   vec ", K[,1]~]), ""))))));}
\q
