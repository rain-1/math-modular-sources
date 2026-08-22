default(parisizemax, 6000000000);
Bn = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn.txt");
ciB = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recB.txt");
QB = vector(36, k, my(i=k-1); sum(j=0,7, ciB[i*8+j+1]*'n^j));
resid(x,n) = { my(sm=0); for(i=0,35, if(n-i>=0, sm += subst(QB[i+1],'n,n)*x[n-i+1])); sm }
NCC = 400;
cn = vector(NCC+1, k, resid(Bn, k-1));
pp = 2^61-1;
cm = vector(NCC+1, k, if(cn[k]==0, Mod(0,pp), Mod(numerator(cn[k]),pp)/Mod(denominator(cn[k]),pp)));
kd(rr,DD) = {
  my(nc = (rr+1)*(DD+1));
  my(nlo = rr+2);
  my(nhi = min(NCC, nlo+nc+40));
  my(nr = nhi-nlo+1);
  if(nr < nc+20, return(-1));
  my(mat = matrix(nr, nc, a, b, cm[nlo+a-1-((b-1)\(DD+1))+1] * Mod(nlo+a-1,pp)^((b-1)%(DD+1))));
  nc - matrank(mat)
}
print("minimal recurrence for c_n = L_A[B]_n :");
{ print("     D: ", vector(13,j,j-1+8));
  for(rr=1,12, print("r=",rr,": ", vector(13,j, kd(rr, j-1+8)))); }
quit
