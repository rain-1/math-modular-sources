/* 04_fricke.gp -- the Atkin-Lehner structure of Phi and the trace identity.
   (a) Phi|_4 W_N = -Phi   (algebraic consequence of u|W_N = 1/(Cu));  numerical check.
   (b) trace identity  Tr^N_{N/p}(Phi) = Phi + p^{1-k/2} Phi|W_p|U_p   (p || N, k=4)
       so if Phi|W_p = eps Phi then  Phi|U_p = (p/eps)(Tr - Phi).
   (c) exact cells: Tr^7_1(Phi_{s7}) = 0 and Tr^10_5(Phi_{s10}) = 0 as q-series.  */
default(parisize, 4000000000);
read("lib.gp");
M = 300; NN = M+10;
{ PH = vector(3); CV = vector(3);
  for(k=1,3, my(S=Setup(k,NN)); PH[k]=S[4]; CV[k]=vector(M,m,polcoeff(S[4],m))); }
Up(f,p,MM) = Ser(concat([0],vector(MM\p,m,polcoeff(f,p*m))),'q);
print("=== (c1) s7:  Tr^7_1(Phi) = Phi - (1/7) Phi|U_7  as a q-series ===");
{ my(T = vector(M\7,m,CV[1][m] - CV[1][7*m]/7));
  print("   coefficients m=1..20 : ",vector(20,m,T[m]));
  print("   all zero for m<=",M\7," ? ",#select(x->x!=0,T)==0); }
print();
print("=== (c2) s10: Tr^10_5(Phi) = Phi - (1/5) Phi|U_5 ===");
{ my(T = vector(M\5,m,CV[2][m] - CV[2][5*m]/5));
  print("   coefficients m=1..20 : ",vector(20,m,T[m]));
  print("   all zero for m<=",M\5," ? ",#select(x->x!=0,T)==0); }
print();
print("=== (c3) s18: Tr^18_9(Phi) = Phi + (1/2) Phi|W_2|U_2 = Phi - (1/2) Phi|U_2  (if Phi|W_2=-Phi) ===");
{ my(T = vector(M\2,m,CV[3][m] - CV[3][2*m]/2));
  print("   coefficients m=1..14 : ",vector(14,m,T[m]));
  print("   is it 0 ? ",#select(x->x!=0,T)==0);
  print("   Tr - 2*Phi mod 4, m=1..14 : ",vector(14,m,lift(Mod(T[m]-2*CV[3][m],4))));
  print("   all Tr = 2 Phi mod 4 for m<=",M\2,"? ",#select(m->(T[m]-2*CV[3][m])%4!=0,vector(M\2,m,m))==0); }
print();
print("=== (c4) s10 at p=2:  Tr^10_5(Phi) = Phi + (1/2)Phi|W_2|U_2;  test Tr = 2 Phi mod 4 with Phi|W_2=+Phi ===");
{ my(T = vector(M\2,m,CV[2][m] + CV[2][2*m]/2));
  print("   Phi + (1/2)Phi|U_2, m=1..14 : ",vector(14,m,T[m]));
  print("   mod 4 vs 2*Phi : ",vector(14,m,lift(Mod(T[m]-2*CV[2][m],4))));
  my(T2 = vector(M\2,m,CV[2][m] - CV[2][2*m]/2));
  print("   Phi - (1/2)Phi|U_2, m=1..14 : ",vector(14,m,T2[m]));
  print("   mod 4 vs 2*Phi : ",vector(14,m,lift(Mod(T2[m]-2*CV[2][m],4)))); }
print();
print("=== (c5) s7 at p=2,5 (p nmid 7): no AL; record Phi|U_p - p Phi over p^3 ===");
{ for(k=1,3, my(row=[]);
   forprime(p=2,13, my(ps=psival(k,p), g=0);
     for(m=1,M\p, g=gcd(g, CV[k][p*m]-ps*p*CV[k][m]));
     row=concat(row,[[p, if(g==0,"inf",valuation(g,p))]]));
   print("  ",NAM[k],": ",row)); }
quit;
