default(realprecision, 90);
OUT="/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/";
{chk(tag, pre) =
  my(A = readvec(Str(OUT,pre,"_A.txt")), NA=#A-1);
  my(R = vector(4, i, my(B=readvec(Str(OUT,pre,"_B",i,".txt"))); B[NA+1]/A[NA+1]*1.0));
  print("\n=== ", tag);
  for(i=1,4, print("   xi_",i," = ", R[i]));
  print("   lindep([xi1..xi4, G])       = ", lindep(concat(R,[Catalan]), 40));
  print("   lindep([xi1..xi4, zeta(2)]) = ", lindep(concat(R,[zeta(2)]), 40));
  print("   lindep([xi1..xi4, G, zeta(2)]) = ", lindep(concat(R,[Catalan, zeta(2)]), 40));}
chk("placement 3 (c=3)", "lvl12_p3");
chk("placement 6 (c=6)", "lvl12_p6");
print("\n-G/27 = ", -Catalan/27, "   G/9 = ", Catalan/9);
quit;
