default(parisize, 6000000000);
default(realprecision, 40);
{mkmat(a,r,D,NF) = matrix(NF, (r+1)*(D+1), i, j,
   my(n = i+r-1, s = (j-1)\(D+1), e = (j-1)%(D+1));
   n^(D-e) * a[n+2-s]);}   /* a[m+1] = a_m */
{doone(fn, tag) =
  my(L = readvec(fn), NA = #L-1);
  print("\n########## ", tag, "   a_0..a_", NA);
  print("  A = ", vector(10,i,L[i]));
  for(r=2,6, for(D=1,4,
     my(NU=(r+1)*(D+1), NF=NU+45);
     if(NF+r+2 > NA, next);
     my(K = matker(mkmat(L,r,D,NF)), kd = matsize(K)[2]);
     if(kd>0,
       print("  >>> terms=",r+1," deg=",D," unknowns=",NU," kerdim=",kd);
       if(kd==1,
         my(v=K[,1], den=1); for(i=1,#v, den=lcm(den,denominator(v[i]))); v*=den;
         my(g=0); for(i=1,#v, g=gcd(g,v[i])); v/=g; if(v[1]<0, v=-v);
         my(P = vector(r+1, s, sum(e=0,D, v[(s-1)*(D+1)+e+1]*nn^(D-e))));
         for(s=1,r+1, print("      P_",s-1,"(n) = ", P[s], "   = ", if(P[s]==0,0,factor(P[s]))));
         my(bad=0); for(n=r, NA-1, if(sum(s=0,r, subst(P[s+1],nn,n)*L[n+2-s])!=0, bad=n;break));
         print("      holds for ALL n=",r,"..",NA-1,"? ", bad==0, if(bad,Str(" first fail n=",bad),""));
         my(lead = vector(r+1,s,polcoeff(P[s],D,nn)));
         my(chp = sum(s=0,r, lead[s+1]*lam^(r-s)));
         print("      char poly ", chp, " = ", factor(chp), "  roots ", polroots(chp));
         return(P))))); 
  0;}
P3 = doone("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl12_p3_A.txt","placement 3 (c=3)");
P6 = doone("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl12_p6_A.txt","placement 6 (c=6)");
quit;
