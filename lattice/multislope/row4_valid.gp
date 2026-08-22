default(parisizemax, 6000000000);
pp = 2^61-1;
NA = 300;
ap = vector(NA+1, i, my(m=i-1); sum(k=0,m, binomial(m,k)^2*binomial(m+k,k)^2));
Am = vector(NA+1, i, Mod(ap[i], pp));
kdim(r, D, nhi) = {
  my(nc = (r+1)*(D+1), nlo = r, nr = nhi-nlo+1);
  if(nr < nc+20, return(-1));
  my(mat = matrix(nr, nc, a, b,
      my(n = nlo + a - 1, i = (b-1)\(D+1), j = (b-1)%(D+1));
      Am[n-i+1] * Mod(n,pp)^j ));
  nc - matrank(mat)
}
print("VALIDATION on Apery zeta(3) numbers (expect first nonzero kernel at r=2,D=3):");
{ for(r=1,4, print("r=",r,": ", vector(7, j, kdim(r, j-1, NA)))); }
quit
