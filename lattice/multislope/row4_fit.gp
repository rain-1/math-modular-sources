default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
NMAX = #An - 1;
pp = 2^61-1;
Am = vector(NMAX+1, i, Mod(An[i], pp));
print("NMAX = ", NMAX);

kdim(r, D, nhi) = {
  my(nc = (r+1)*(D+1), nlo = r, nr = nhi-nlo+1);
  if(nr < nc+5, return(-1));
  my(mat = matrix(nr, nc, a, b,
      my(n = nlo + a - 1, i = (b-1)\(D+1), j = (b-1)%(D+1));
      Am[n-i+1] * Mod(n,pp)^j ));
  nc - matrank(mat)
}

{ print("kernel dimension table: rows r=1..10, cols D=0..12  (-1 = not enough rows)");
  print("      D:  ", vector(13,j,j-1));
  for(r=1,10,
    print("r=",r,": ", vector(13, j, kdim(r, j-1, min(NMAX, 400))))
  );
}
quit
