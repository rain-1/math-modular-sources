default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
Bn = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn.txt");
NMAX = #An - 1;
pp = 2^61-1;
Am = vector(NMAX+1, i, Mod(An[i],pp));
Bm = vector(NMAX+1, i, Mod(numerator(Bn[i]),pp)/Mod(denominator(Bn[i]),pp));
kdimj(r, D, extra) = {
  my(nc = (r+1)*(D+1));
  my(nlo = r);
  my(nhi = min(NMAX, nlo+(nc\2)+extra));
  my(nr = nhi-nlo+1);
  if(2*nr < nc+20, return(-1));
  my(m1 = matrix(nr, nc, a, b, Am[nlo+a-1-((b-1)\(D+1))+1] * Mod(nlo+a-1,pp)^((b-1)%(D+1))));
  my(m2 = matrix(nr, nc, a, b, Bm[nlo+a-1-((b-1)\(D+1))+1] * Mod(nlo+a-1,pp)^((b-1)%(D+1))));
  nc - matrank(matconcat([m1;m2]))
}
print("JOINT kernel dims (operator annihilating BOTH A_n and B_n):");
{ print("     D: ", vector(8,j,j+5));
  for(r=28,50, print("r=",r,": ", vector(8,j, kdimj(r, j+5, 80)))); }
quit
