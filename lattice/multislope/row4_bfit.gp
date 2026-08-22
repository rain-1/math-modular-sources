default(parisizemax, 6000000000);
Bn = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn.txt");
NMAX = #Bn - 1;
pp = 2^61-1;
Bm = vector(NMAX+1, i, Mod(numerator(Bn[i]),pp)/Mod(denominator(Bn[i]),pp));
print("NMAX = ", NMAX);
kdim(r, D, extra) = {
  my(nc = (r+1)*(D+1));
  my(nlo = r);
  my(nhi = min(NMAX, nlo+nc+extra));
  my(nr = nhi-nlo+1);
  if(nr < nc+20, return(-1));
  my(mat = matrix(nr, nc, a, b, Bm[nlo+a-1-((b-1)\(D+1))+1] * Mod(nlo+a-1,pp)^((b-1)%(D+1))));
  nc - matrank(mat)
}
{ print("        D: ", vector(10,j,j+5));
  for(r=22,50, print("r=",r,": ", vector(10, j, kdim(r, j+5, 60))));
}
quit
