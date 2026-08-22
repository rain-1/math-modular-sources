default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An7.txt");
NMAX = #An - 1;
pp = 2^61-1;
Am = vector(NMAX+1, i, Mod(An[i], pp));
kdim(r, D, extra) = {
  my(nc = (r+1)*(D+1), nlo = r, nhi = min(NMAX, nlo+nc+extra), nr = nhi-nlo+1);
  if(nr < nc+20, return(-1));
  my(mat = matrix(nr, nc, a, b,
      my(n = nlo + a - 1, i = (b-1)\(D+1), j = (b-1)%(D+1));
      Am[n-i+1] * Mod(n,pp)^j ));
  nc - matrank(mat)
}
{ print("rows r (shift order) = 14..40, cols D (poly degree) = 5..14");
  print("        D: ", vector(10,j,j+4));
  for(r=14,40,
    print("r=",if(r<10," ",""),r,": ", vector(10, j, kdim(r, j+4, 40)));
  );
}
quit
