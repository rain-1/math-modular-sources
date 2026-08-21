read("lattice/euler_criterion/rows.gp");
read("lattice/theoremF_hyp/vlaw.gp");
MM = 3000;
{
  my(rs);
  rs = [
   ["A (7,2,-8)",   rowR2(7,2,-8,MM)[1],   2],
   ["B (9,3,27)",   rowR2(9,3,27,MM)[1],   3],
   ["C (10,3,9)",   rowR2(10,3,9,MM)[1],   3],
   ["E (12,4,32)",  rowR2(12,4,32,MM)[1],  2],
   ["F (17,6,72)",  rowR2(17,6,72,MM)[1],  2],
   ["F (17,6,72)",  rowR2(17,6,72,MM)[1],  3],
   ["alpha Domb (10,4,64)", rowR3(10,4,64,MM)[1], 2],
   ["delta (7,3,81)",  rowR3(7,3,81,MM)[1],  3],
   ["eps (12,4,16)",   rowR3(12,4,16,MM)[1], 2],
   ["zeta (9,3,-27)",  rowR3(9,3,-27,MM)[1], 3],
   ["eta (11,5,125)",  rowR3(11,5,125,MM)[1],5],
   ["s18 Cooper",      rowS18(MM)[1],        3],
   ["cusp L12 wt3",    rowCusp(MM)[1],       2]
  ];
  for(i=1,#rs, vlaw(rs[i][2], rs[i][3], rs[i][1], MM); super(rs[i][2], rs[i][3], rs[i][1], MM); print());
}
quit;
