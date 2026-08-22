read("05_jtest.gp");
tst(nm,M,j1,j2,A,B,C) = {
  my(r = jtest(M,j1,j2,A,B,C));
  if(r==0, printf("%-26s (%d;%d,%d) A=%d B=%d C=%d  ->  monodromy NOT in PSL2(Z)\n",nm,M,j1,j2,A,B,C),
           printf("%-26s (%d;%d,%d) A=%d B=%d C=%d  ->  h=%d  deg J=%d   J=%s / %s\n",nm,M,j1,j2,A,B,C,r[1],r[2],r[3],r[4]));
}
tst("Zagier A (7,2,-8)", 1,0,0, 7,2,-8);
tst("Zagier B (9,3,27)", 1,0,0, 9,3,27);
tst("Zagier C (10,3,9)", 1,0,0, 10,3,9);
tst("Zagier D (11,3,-1)",1,0,0, 11,3,-1);
tst("Zagier E (12,4,32)",1,0,0, 12,4,32);
tst("Zagier F (17,6,72)",1,0,0, 17,6,72);
tst("root Apery",        2,1,1, 136,10,4);
tst("root T",            2,1,1, 24,2,4);
tst("root Domb",         2,1,1, 20,2,16);
tst("root AZ(9,3,-27)",  2,1,1, 72,6,-108);
tst("root AZ(11,5,125)", 2,1,1, 88,10,500);
tst("root AZ(7,3,81)",   2,1,1, 56,6,324);
tst("root Cooper s7",    3,1,2, 26,2,-3);
quit;
