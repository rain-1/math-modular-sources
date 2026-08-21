default(realprecision, 40);
read("lattice/theoremF_hyp/descent.gp");
/* number of cusps with denominator c on X_0(N) is eulerphi(gcd(c,N/c)) */
degt2(N,L) = { my(s=0); fordiv(N,c, my(o=ligoz(N,L,c)); if(o<0, s -= o*eulerphi(gcd(c,N/c)))); s; }
rep2(name,N,Lt,LF,wt,p) = { print("--- ",name,"  N=",N," deg t (with cusp multiplicity) = ", degt2(N,Lt));
  report(name,N,Lt,LF,wt,p); }
rep2("B twin level 9", 9, [[9,3],[1,-3]], [[1,3],[3,-1]], 1, 3);
rep2("eta twin level 5", 5, [[5,6],[1,-6]], [[1,5],[5,-1]], 2, 5);
rep2("delta twin level 6", 6, [[3,4],[6,4],[1,-4],[2,-4]], [[1,3],[2,3],[3,-1],[6,-1]], 2, 3);
rep2("cusp row level 12", 12, [[1,6],[3,6],[4,6],[12,6],[2,-6],[6,-6]], [[2,5],[6,5],[1,-2],[3,-2],[4,-2],[12,-2]], 1, 2);
rep2("cusp row level 6 twin", 6, [[2,6],[6,6],[1,-6],[3,-6]], [[1,2],[3,2],[2,-1],[6,-1]], 1, 2);
/* recompute all canonical rows with cusp multiplicity */
print("### deg t with cusp multiplicity, canonical presentations");
print("A  ", degt2(6,[[1,3],[2,-3],[3,-9],[6,9]]));
print("C  ", degt2(6,[[1,4],[2,-8],[3,-4],[6,8]]));
print("E  ", degt2(8,[[1,4],[2,-10],[4,2],[8,4]]));
print("F  ", degt2(12,[[1,5],[2,-14],[3,1],[4,5],[6,2],[12,1]]));
print("B  ", degt2(36,[[1,3],[2,-9],[4,3],[9,-3],[18,9],[36,-3]]));
print("al ", degt2(12,[[1,6],[2,-12],[3,6],[4,6],[6,-12],[12,6]]));
print("ga ", degt2(6,[[1,12],[2,-12],[3,-12],[6,12]]));
print("de ", degt2(12,[[1,4],[2,-16],[3,-4],[4,4],[6,16],[12,-4]]));
print("ep ", degt2(8,[[1,8],[2,-8],[4,-8],[8,8]]));
print("ze ", degt2(9,[[1,6],[3,-12],[9,6]]));
print("et ", degt2(20,[[1,6],[2,-18],[4,6],[5,-6],[10,18],[20,-6]]));
quit;
