read("05_jtest.gp");
tst(nm,M,j1,j2,A,B,C) = {
  my(r = jtest(M,j1,j2,A,B,C));
  if(r==0, printf("%-24s -> NOT in PSL2(Z)\n",nm),
           printf("%-24s -> h=%d  degJ=%d  gamma=%s\n",nm,r[1],r[2],r[3]));
}
tst("Zagier A (7,2,-8)", 1,0,0, 7,2,-8);
tst("Zagier C (10,3,9)", 1,0,0, 10,3,9);
tst("Zagier D (11,3,-1)",1,0,0, 11,3,-1);
tst("Zagier E (12,4,32)",1,0,0, 12,4,32);
tst("Zagier F (17,6,72)",1,0,0, 17,6,72);
quit;
