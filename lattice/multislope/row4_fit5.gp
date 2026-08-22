default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
NMAX = #An - 1;
pp = 2^61-1;
Am = vector(NMAX+1, i, Mod(An[i], pp));
print("NMAX = ", NMAX);
kdim(r, D, extra) = {
  my(nc = (r+1)*(D+1), nlo = r, nhi = min(NMAX, nlo+nc+extra), nr = nhi-nlo+1);
  if(nr < nc+20, return(-1));
  my(mat = matrix(nr, nc, a, b,
      my(n = nlo + a - 1, i = (b-1)\(D+1), j = (b-1)%(D+1));
      Am[n-i+1] * Mod(n,pp)^j ));
  nc - matrank(mat)
}
{ print("Is there ANY operator with shift order r and poly degree <= 35?");
  for(r=10,24, print("r=",r,"  D<=35 : kdim = ", kdim(r,35,60)));
  print("r=23  D<=40 : kdim = ", kdim(23,40,60));
  print("r=24  D<=11 : kdim = ", kdim(24,11,300));
  print("r=24  D<=12 : kdim = ", kdim(24,12,300));
  print("r=35  D<=7  : kdim = ", kdim(35,7,600));
  print("r=34  D<=7  : kdim = ", kdim(34,7,600));
}
quit
