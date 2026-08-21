read("lattice/theoremF_hyp/twins.gp");
read("lattice/euler_criterion/rows.gp");
N0=24;
/* row B: (9,3,27) */
aB = rowR2(9,3,27,N0)[1];
aBflip = vector(N0+1, i, (-1)^(i-1)*aB[i]);
t9 = etaq([[9,3],[1,-3]]);  F9 = etaq([[1,3],[3,-1]]);
rowcheck("B level 9 (t=eta9^3/eta1^3, F=eta1^3/eta3)", t9, F9, aBflip, 20);
print();
/* row eta: (11,5,125) */
aE = rowR3(11,5,125,N0)[1];
aEflip = vector(N0+1, i, (-1)^(i-1)*aE[i]);
t5 = etaq([[5,6],[1,-6]]);  F5 = etaq([[1,5],[5,-1]]);
rowcheck("eta level 5 (t=eta5^6/eta1^6, F=eta1^5/eta5)", t5, F5, aEflip, 20);
print();
/* row delta: (7,3,81) */
aD = rowR3(7,3,81,N0)[1];
aDflip = vector(N0+1, i, (-1)^(i-1)*aD[i]);
t6 = etaq([[3,4],[6,4],[1,-4],[2,-4]]);  F6 = etaq([[1,3],[2,3],[3,-1],[6,-1]]);
rowcheck("delta level 6 (t=e3^4e6^4/e1^4e2^4, F=e1^3e2^3/e3e6)", t6, F6, aDflip, 20);
print();
/* sources */
print("Phi_B' = F9*Dq(t9) : ", Vec(Pol(F9*Dq(t9)+O(q^13))));
print("Phi_eta' = F5*Dq(t5) : ", Vec(Pol(F5*Dq(t5)+O(q^13))));
print("Phi_delta' = F6*Dq(t6) : ", Vec(Pol(F6*Dq(t6)+O(q^13))));
quit;
