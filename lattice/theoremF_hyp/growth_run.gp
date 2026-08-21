read("lattice/euler_criterion/rows.gp");
Nlim = 3000;
NN=3000;
read("lattice/theoremF_hyp/growth.gp");
{
  my(rs);
  rs = [
   ["A (7,2,-8)",   rowR2(7,2,-8,NN)[1],   2],
   ["B (9,3,27)",   rowR2(9,3,27,NN)[1],   3],
   ["C (10,3,9)",   rowR2(10,3,9,NN)[1],   3],
   ["D (11,3,-1)",  rowR2(11,3,-1,NN)[1],  2],
   ["E (12,4,32)",  rowR2(12,4,32,NN)[1],  2],
   ["F (17,6,72)",  rowR2(17,6,72,NN)[1],  2],
   ["F (17,6,72)",  rowR2(17,6,72,NN)[1],  3],
   ["alpha Domb (10,4,64)", rowR3(10,4,64,NN)[1], 2],
   ["gamma Apery (17,5,1)", rowR3(17,5,1,NN)[1],  2],
   ["delta (7,3,81)",  rowR3(7,3,81,NN)[1],  3],
   ["eps (12,4,16)",   rowR3(12,4,16,NN)[1], 2],
   ["zeta (9,3,-27)",  rowR3(9,3,-27,NN)[1], 3],
   ["eta (11,5,125)",  rowR3(11,5,125,NN)[1],5],
   ["s18 Cooper",      rowS18(NN)[1],        3],
   ["cusp L12 wt3",    rowCusp(NN)[1],       2]
  ];
  for(i=1,#rs, profile(rs[i][2], rs[i][3], rs[i][1], NN); print());
}
quit;
