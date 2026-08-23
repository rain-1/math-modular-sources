/* driver: the single-congruence lattice {c : c.(Y,U) = 0 mod S_n} of 06_threshold.tex */
\p 3000
default(parisize, 2000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW = rdrows(concat(DIR,"rows_all.txt"));
print("tag,n,logidx,halflogcov,logl1,logl2,skew,incone,logcone,rho,rho2,gcdflag,cflen,cfmax,cfidx,cfkind,cfpq");
for(n=4,120, analS(n, GG, RW));
\q
