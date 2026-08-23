/* driver: full structural table, n = 4..NHI, k in KLIST.
   cat lattice/positivity/rows_pos.gp lattice/p2_structure/{p2core.gp,p2run.gp} \
       lattice/p2_structure/run_main.gp > run.gp && gp -q run.gp                */
\p 3000
default(parisize, 2000000000);
GG = Catalan;
NHI = 120;
KLIST = [22.4, 23.0, 23.9];
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW = rdrows(concat(DIR,"rows_all.txt"));
DUMP = concat(DIR,"vectors_n120.txt");
print("k,n,logidx,halflogcov,logl1,logl2,skew,incone,logcone,rho,rho2,logq,jhi,mrg,ok,cflen,cfmax,cfidx,cfkind,cfpq,logh11,h22,logsinth,pqnext,cvkind,cvidx,cvt");
for(n=4,NHI, anal(n, KLIST, GG, RW, DUMP));
\q
