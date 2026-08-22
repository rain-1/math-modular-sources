read("05_jtest.gp");
tst(nm,M,j1,j2,A,B,C) = {
  my(r = jtest(M,j1,j2,A,B,C));
  if(r==0, printf("%-24s -> monodromy NOT conjugate into PSL2(Z)\n",nm),
           printf("%-24s -> h=%d  deg J=%d  gamma=%s\n",nm,r[1],r[2],r[3]));
}
tst("root AZ(9,3,-27)",  2,1,1, 72,6,-108);
tst("root AZ(11,5,125)", 2,1,1, 88,10,500);
tst("root Domb",         2,1,1, 20,2,16);
tst("root Apery",        2,1,1, 136,10,4);
tst("root T",            2,1,1, 24,2,4);
tst("root AZ(7,3,81)",   2,1,1, 56,6,324);
tst("root Cooper s7",    3,1,2, 26,2,-3);
tst("Zagier B (9,3,27)", 1,0,0, 9,3,27);
quit;
