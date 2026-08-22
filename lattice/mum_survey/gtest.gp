read("gammap.gp");
read("../euler_criterion/lp.gp");
{
for(pi=1,3,
  my(p=[2,3,5][pi], PR=60, K=14);
  my(c = gpTaylor(p,K,PR));
  print("p=",p);
  for(k=1,5, print("   c_",k," = ", c[k]));
  my(z2=Lp(p,triv,2,40), z3=Lp(p,triv,3,40));
  print("   zeta_p(2) = ",z2);
  print("   zeta_p(3) = ",z3);
  my(L2 = 2*c[2], L3 = 6*c[3]);
  print("   L2=2c2 = ",L2, "   L2/z2 = ", L2/z2);
  print("   L3=6c3 = ",L3, "   L3/z3 = ", L3/z3);
);
}
quit
