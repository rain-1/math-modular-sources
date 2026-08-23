/* driver: Task 2 -- the rationality-kernel vector.
   cat lattice/positivity/rows_pos.gp lattice/p2_structure/{p2core.gp,kernel.gp} \
       lattice/p2_structure/run_kernel.gp > run.gp && gp -q run.gp            */
\p 3000
default(parisize, 2000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW = rdrows(concat(DIR,"rows_all.txt"));
print("== kernel identity for rational surrogates ==");
print("e,k,n,zerof,equal,inlat,logcoordZ,logpred,t,g,logkerlen,logpredlen");
{for(e=2,6, my(GS=bestappr(GG,10^e)); for(n=6,24, kerchk(n,[22.4,23.9],GS,e,RW)));}
print("== critical denominator in the true-G metric ==");
print("k,n,logbcrit,halflogcov,logl1,logl2,logMODoveridx");
for(n=4,120, bcrit(n,[22.4,23.0,23.9],GG,RW));
\q
