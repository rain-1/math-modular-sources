default(parisize, 6000000000);
{mkmat(a,r,D,NF) = matrix(NF, (r+1)*(D+1), i, j,
   my(n = i+r-1, s = (j-1)\(D+1), e = (j-1)%(D+1)); n^(D-e) * a[n+2-s]);}
{ck(fn,tag) = my(L=readvec(fn), NA=#L-1);
  print("\n### ",tag," (a_0..a_",NA,")");
  for(r=4,4, for(D=1,8, my(NU=(r+1)*(D+1), NF=NU+50);
     if(NF+r+2>NA, print("   r=",r," D=",D,": not enough terms"); next);
     print("   5-term, deg ",D,": kerdim = ", matsize(matker(mkmat(L,r,D,NF)))[2])));}
ck("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl12_p3_A.txt","placement 3");
ck("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl12_p6_A.txt","placement 6");
quit;
