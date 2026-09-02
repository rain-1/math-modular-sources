\\ 33_halfint.gp -- Task (2): independent construction of PZ's weight-5/2 inputs
\\ G0,G1,G2 at level 4, Lemma 1's F4A,F4B, validation of a(m^2)=m*beta(m).
default(parisize, 6000000000);
NP = 2600;
BP = NP + 40;
TH = 1 + 2*sum(n=1, sqrtint(BP), 'q^(n^2)) + O('q^BP);
{ E24 = sum(n=1, BP-1, if(n%2, sigma(n,1), 0)*'q^n) + O('q^BP); }
{ E4q(t,P) = 1 + 240*sum(n=1, (P-1)\t, sigma(n,3)*'q^(t*n)) + O('q^P); }
{ E6q(t,P) = 1 - 504*sum(n=1, (P-1)\t, sigma(n,5)*'q^(t*n)) + O('q^P); }
E44 = E4q(4,BP);
E64 = E6q(4,BP);
DEL4 = (E44^3 - E64^2)/1728;
G0 = TH*(TH^4 - 20*E24);
G1 = TH*E44^2*E64/DEL4;
G2 = G0*E44^3/DEL4;
print("=== check against PZ's printed expansions (section 4) ===");
{ show(nam,g,lst) = print("  ",nam,": ", vector(#lst,i, Str("a(",lst[i],")=",polcoeff(g,lst[i])))); }
show("G0", G0, [0,1,4,5,8,9,16,25,36,49]);
show("G1", G1, [-4,-3,0,4,9,16,25,36,49]);
show("G2", G2, [-4,-3,0,1,4,9,16,25,36,49]);
print("  PZ g0 : 1, -10, -70, -48, -120, -250, -550, -1210, -1750, -3370");
print("  PZ g1 : 1, 2, 2, -196884, -85975040, -86169224844, -51186246451200, -35015148280961780, -21434928162930081792");
print("  PZ g2 : 1, -10, 674, -7488, 144684, -224574272, -42882054732, -63793268216640, -31501841125150388, -22385069000981561664");
F4A = 7/8*G0 + 1/768*G1 - 1/768*G2;
F4B = 19/18*G0 - 5/648*G1 - 1/648*G2;
print("");
print("=== Lemma 1 ===");
print("  F4A low coefficients a(-4..4) : ", vector(9,i, Str(i-5,":",polcoeff(F4A,i-5))));
print("  F4B low coefficients a(-4..4) : ", vector(9,i, Str(i-5,":",polcoeff(F4B,i-5))));
print("  64*F4A integral (n<=400)? ", if(sum(n=-4,400, denominator(64*polcoeff(F4A,n))!=1)==0,"YES","NO"));
print("  108*F4B integral (n<=400)? ", if(sum(n=-4,400, denominator(108*polcoeff(F4B,n))!=1)==0,"YES","NO"));
{ plusok(f) = my(b=0); for(n=-4,600, if(!(n%4==0||n%4==1), if(polcoeff(f,n)!=0, b=n; break))); b; }
print("  plus space (a(n)=0 unless n=0,1 mod 4, n<=600): F4A ", if(plusok(F4A)==0,"YES","NO"), "  F4B ", if(plusok(F4B)==0,"YES","NO"));
NL = 60;
E4s = E4q(1,NL+10);
E6s = E6q(1,NL+10);
DELs = (E4s^3 - E6s^2)/1728;
LIFTa = DELs/E4s^2;
LIFTb = E4s*DELs/E6s^2;
{ betaof(F,M) = my(A=vector(M,m,polcoeff(F,m)), cp, b);
  cp = vector(M,n,A[n]/n); b=vector(M);
  for(n=1,M, my(s=0); fordiv(n,d, s+=moebius(d)*cp[n/d]); b[n]=s); b; }
B4a = betaof(LIFTa,NL);
B4b = betaof(LIFTb,NL);
print("");
print("=== DICTIONARY CHECK:  a(m^2) = m * beta(m),  beta = (A/n) * mu ===");
{ for(m=1,20,
   my(x=polcoeff(F4A,m^2), y=m*B4a[m], u=polcoeff(F4B,m^2), v=m*B4b[m]);
   print("  m=",m,"  F4a a(m^2)=",x,"  ",if(x==y,"OK","MISMATCH"),"   F4b a(m^2)=",u,"  ",if(u==v,"OK","MISMATCH")));
}
print("  F4a all m<=20 agree: ", if(sum(m=1,20, polcoeff(F4A,m^2)!=m*B4a[m])==0,"YES","NO"));
print("  F4b all m<=20 agree: ", if(sum(m=1,20, polcoeff(F4B,m^2)!=m*B4b[m])==0,"YES","NO"));
print("");
print("=== lift formula A(n) = sum_{d|n} d a(n^2/d^2) reproduces the level-one forms ===");
{ for(j=1,2, my(f=[F4A,F4B][j], F=[LIFTa,LIFTb][j], bad=0);
   for(n=1,50, my(s=0); fordiv(n,d, s += d*polcoeff(f,(n/d)^2)); if(s!=polcoeff(F,n), bad=n; break));
   print("  ", ["F4a","F4b"][j], " : ", if(bad==0,"lift matches for n<=50",Str("FAIL at n=",bad))));
}
write("33_f4a.txt", vector(NP+5, i, polcoeff(F4A, i-5)));
write("33_f4b.txt", vector(NP+5, i, polcoeff(F4B, i-5)));
quit;
