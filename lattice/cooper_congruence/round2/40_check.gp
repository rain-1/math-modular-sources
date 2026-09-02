read("40_core.gp");
N=12;
print("A_n (n=0..12):");
{ for(k=1,3, print("  ",NAM[k],": ",AVEC(k,N))); }
print();
print("a_j (j=0..9):");
{ for(k=1,3, print("  ",NAM[k],": ",vector(10,j,AJ(k,N)[j]))); }
print();
print("check A_n = [x^n]F against q-expansion:");
{ for(k=1,3, my(S=Setup(k,40), F=S[2], x=S[3], A=AVEC(k,20), R, ok);
   /* invert x(q) numerically as series: express F in terms of x */
   R = subst(Ser(A,'t), 't, x + O('q^25));
   ok = 1;
   for(m=0,20, if(polcoeff(R-F,m)!=0, ok=0));
   print("  ",NAM[k],": ",if(ok,"MATCH","MISMATCH"))); }
quit;
