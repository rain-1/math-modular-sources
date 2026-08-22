read("../mum_survey/ops.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
P = OP[4];
ck(k) = sum(i=1,#P, polcoef(P[i],k,'X)*z^(i-1));
S(k,j) = sum(r=0,j, (-1)^(j-r)*binomial(j,r)*r^k)/j!;
{ Q = vector(5, jj, my(j=jj-1); z^j*sum(k=0,4, if(j==0, if(k==0,1,0), S(k,j))*ck(k))); }
w = Mod('w, 'w^2-17);
{ forstep(sg=-1,1,2,
  my(z0 = (349 + sg*85*w)/131072);
  my(qs = vector(5,jj, subst(Q[jj], z, z0+t)));
  my(ords = vector(5,jj, if(qs[jj]==0, 10^9, valuation(qs[jj], t))));
  my(mu = 10^9); for(jj=1,5, if(ords[jj]<10^8, mu=min(mu, ords[jj]-(jj-1))));
  my(ind = 0);
  for(jj=1,5, if(ords[jj]-(jj-1)==mu,
     my(lc=polcoef(qs[jj],ords[jj],t), j=jj-1); ind += lc*prod(r=0,j-1,rho-r)));
  print("z0 = (349 ", if(sg>0,"+","-"), " 85 sqrt17)/131072   mu=",mu);
  print("   indicial = ", ind);
  
  print("   exponents (numeric) = ", polroots(subst(lift(ind), 'w, sqrt(17.))));
  ); }
quit
