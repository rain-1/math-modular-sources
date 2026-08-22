default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
NMAX = #An - 1;
pp = 2^61-1;
Am = vector(NMAX+1, i, Mod(An[i], pp));

kdim(r, D, nhi) = {
  my(nc = (r+1)*(D+1), nlo = r, nr = nhi-nlo+1);
  if(nr < nc+20, return(-1));
  my(mat = matrix(nr, nc, a, b,
      my(n = nlo + a - 1, i = (b-1)\(D+1), j = (b-1)%(D+1));
      Am[n-i+1] * Mod(n,pp)^j ));
  nc - matrank(mat)
}
{ my(Dmax=25);
  print("      D:  ", vector(Dmax+1,j,j-1));
  for(r=1,20,
    print("r=",r,": ", vector(Dmax+1, j, kdim(r, j-1, NMAX)))
  );
}
quit
