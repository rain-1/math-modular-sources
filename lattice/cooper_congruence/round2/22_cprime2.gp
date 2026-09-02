\\ 22_cprime2.gp -- consequence of the mod-2 law: c'(m) mod 2 is multiplicative and
\\ c'(m) is odd  <=>  the largest divisor of m coprime to 2P (resp. to 2 and P) is a square.
read("lib.gp");
CP = [read("20b_cp_s7.txt"), read("20b_cp_s10.txt"), read("20b_cp_s18.txt")];
G  = [read("20b_gamma_s7.txt"), read("20b_gamma_s10.txt"), read("20b_gamma_s18.txt")];
N  = #CP[1];
print("N = ", N);
PSP = [7,5,3];
print();
print("=== c'(m) mod 2 : predicted odd  <=>  the 2P-free part of m is a perfect square ===");
{
for(k=1,3,
  my(P=PSP[k], bad);
  bad = select(m->my(t=m); while(t%2==0, t=t/2); while(t%P==0, t=t/P);
                (CP[k][m]%2==1) != issquare(t),
              vector(N,i,i));
  print("  ",NAM[k]," (P=",P,"): ", if(#bad==0,"HOLDS",concat("FAILS first at m=",bad[1])));
);
}
print();
print("=== c'(m) mod 2 multiplicative on coprime pairs? ===");
{
for(k=1,3,
  my(f=0,t=0);
  for(a=2,80, for(b=2,N\a, if(gcd(a,b)==1, t++; if((CP[k][a*b]-CP[k][a]*CP[k][b])%2!=0, f++))));
  print("  ",NAM[k],": failures ",f,"/",t);
);
}
print();
print("=== unified form of the mod-2 law: gamma(n) odd <=> 2^{1+v_2(N)} nmid n and P nmid n ===");
{
for(k=1,3,
  my(P=PSP[k], a=1+valuation(LEV[k],2), bad);
  bad = select(n->((G[k][n]%2==1) != (n%(2^a)!=0 && n%P!=0)), vector(N,i,i));
  print("  ",NAM[k],": N=",LEV[k],", 2^",a,"=",2^a,", P=",P," : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
);
}
quit;
