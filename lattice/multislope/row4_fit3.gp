default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An7.txt");
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
{ for(D=6,10,
    print("D=",D,": ", vector(11, k, my(r=20+5*(k-1)); kdim(r,D,60)));
  );
}
quit
