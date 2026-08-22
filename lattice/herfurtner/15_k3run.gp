read("/home/ubuntu/code/math-modular-sources/lattice/herfurtner/15_k3jtest.gp");

print("K3 window: dmax=", DMAX, " hmax=", HMAX, " terms=", NTERM, " #gamma=", 4*#GBASE, " p=", KPRIME);
print("label | M j1 j2 A B C | verdict h degJ gamma used extra | ms");
tst(nm, mm, j1, j2, A, B, C) =
{ my(r, ms);
  gettime(); r = k3test(mm, j1, j2, A, B, C); ms = gettime();
  if(r == 0,
     printf("%-42s | %d %d %d %d %d %d | NO  -  -  -  -  - | %d\n", nm, mm, j1, j2, A, B, C, ms),
     printf("%-42s | %d %d %d %d %d %d | YES %d %d %s %d %d | %d\n", nm, mm, j1, j2, A, B, C, r[1], r[2], r[3], r[4], r[5], ms));
}

tst("CONTROL+ Zagier B", 1, 0, 0, 9, 3, 27);
tst("CONTROL+ Zagier D", 1, 0, 0, 11, 3, -1);
tst("CONTROL+ sqrt(AZ(11,5,125))        [root row]", 2, 1, 1, 88, 10, 500);
tst("CONTROL+ sqrt(AZ(9,3,-27))         [root row]", 2, 1, 1, 72, 6, -108);
tst("CONTROL+ NEW row on I1 I7 II II", 3, 1, 1, 117, 21, 441);
tst("CONTROL+ NEW row on I3 III III III", 4, 1, 3, 72, 6, 108);
tst("sec6.1-NEG sqrt(Domb)                [root row]", 2, 1, 1, 20, 2, 16);
tst("sec6.1-NEG sqrt(T)                   [root row]", 2, 1, 1, 24, 2, 4);
tst("sec6.1-NEG sqrt(AZ(7,3,81))          [root row]", 2, 1, 1, 56, 6, 324);
tst("sec6.1-NEG Beukers = sqrt(Apery)     [root row]", 2, 1, 1, 136, 10, 4);
tst("sec6.1-NEG A=0 row, class (2;1,3) -", 2, 1, 3, 0, -4, -64);
tst("sec6.1-NEG A=0 row, class (2;1,3) +", 2, 1, 3, 0, 4, -64);
tst("sec6.1-NEG new row (3;1,1), not a surface", 3, 1, 1, 180, 24, -72);
tst("sec6.1-NEG sqrt(s7) Cooper           [root row]", 3, 1, 2, 26, 2, -3);
tst("sec6.1-NEG new row (3;1,2), not a surface", 3, 1, 2, 40, 4, 48);
tst("sec6.1-NEG A=0 row, class (3;2,4) -", 3, 2, 4, 0, -3, -81);
tst("sec6.1-NEG A=0 row, class (3;2,4) +", 3, 2, 4, 0, 3, -81);
tst("sec6.1-NEG new row (4;1,3), not a surface", 4, 1, 3, 28, 2, -8);
tst("sec6.1-NEG new row (4;1,3), not a surface", 4, 1, 3, 48, 4, 32);
tst("sec6.1-NEG new row (4;1,3), not a surface", 4, 1, 3, 68, 6, 72);
tst("sec6.1-NEG new row (4;1,3), not a surface", 4, 1, 3, 80, 6, 36);
tst("sec6.1-NEG level-5 Fricke row", 4, 1, 3, 88, 6, -4);
tst("sec6.1-NEG A=0 row, class (4;3,5) -", 4, 3, 5, 0, -4, -256);
tst("sec6.1-NEG A=0 row, class (4;3,5) +", 4, 3, 5, 0, 4, -256);
tst("CONTROL- sqrt(s10)   [not Kodaira: delta_inf = 1/4]", 8, 3, 5, 24, 2, -4);
tst("CONTROL- sqrt(s18)   [not Kodaira: delta_inf = 1/4]", 8, 3, 5, 56, 6, 12);
tst("CONTROL- NCS 4.4 arctan row, alpha=18  [rho = 7/6]", 3, 2, 5, 18, -6, -3);
quit;
