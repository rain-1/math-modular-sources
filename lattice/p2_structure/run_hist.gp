/* driver: Gauss-Kuzmin histogram of the partial quotients of h12/h11 */
\p 3000
default(parisize, 2000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data/";
RW = rdrows(concat(DIR,"rows_all.txt"));
print("k,n,npq,a1,a2,a3,a4,a5,a6,a7,a8,a9plus");
for(n=4,120, cfhist(n, [22.4,23.0,23.9], GG, RW));
\q
