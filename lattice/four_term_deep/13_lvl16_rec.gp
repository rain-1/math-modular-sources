default(parisize, 3000000000);
L  = readvec("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_A.txt");
NA = #L - 1; a(n) = L[n+1];
{mkmat(r,D,NF) = matrix(NF, (r+1)*(D+1), i, j,
   my(n = i+r-1, s = (j-1)\(D+1), e = (j-1)%(D+1));
   n^(D-e) * a(n+1-s));}
{show(r,D) =
  my(NF = (r+1)*(D+1)+45, K = matker(mkmat(r,D,NF)));
  print("\n===== terms=",r+1," deg=",D,"  kerdim=",matsize(K)[2]);
  if(matsize(K)[2]!=1, return(0));
  my(v = K[,1], den = 1); for(i=1,#v, den = lcm(den, denominator(v[i])));
  v = v*den; my(g = 0); for(i=1,#v, g = gcd(g, v[i])); v = v/g;
  if(v[1]<0, v = -v);
  my(P = vector(r+1, s, sum(e=0,D, v[(s-1)*(D+1)+e+1]*nn^(D-e))));
  for(s=1,r+1, print("  P_",s-1,"(n) [coef of a_{n+1-",s-1,"}] = ", P[s], "   factored: ", if(P[s]==0,0,factor(P[s]))));
  /* verify over all n */
  my(bad=0);
  for(n=r, NA-1,
    my(t = sum(s=0,r, subst(P[s+1],nn,n)*a(n+1-s)));
    if(t!=0, bad=n; break));
  print("  holds for ALL n=",r,"..",NA-1,"? ", bad==0, if(bad,Str("  first failure n=",bad),""));
  /* characteristic polynomial: leading coefficients in n */
  my(lead = vector(r+1, s, polcoeff(P[s], D, nn)));
  my(chp = sum(s=0,r, lead[s+1]*lam^(r-s)));
  print("  leading-coeff vector (in n^",D,"): ", lead);
  print("  char poly (in lam, from leading coeffs): ", chp, "  -> factor ", factor(chp));
  print("  roots lam (=1/x_sing): ", polroots(chp));
  print("  reciprocals x = 1/lam: ", vector(poldegree(chp), i, 1/polroots(chp)[i]));
  P;
}
default(realprecision, 40);
R62 = show(5,2);
R53 = show(4,3);
