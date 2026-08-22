default(parisizemax, 8000000000);
outdir = "/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/out/";
poch(a,m) = prod(k=0,m-1, a+k);
qc6(m) = { 135/2 * (-1024)^m * poch(1/2,m)*poch(1/12,m)^2*poch(5/12,m)^2*poch(7/12,m)^2*poch(11/12,m)^2 / (poch(1,m)^3*poch(1/3,m)*poch(5/3,m)*poch(1/6,m)^2*poch(5/6,m)^2); }
qc12(m) = { 1152 * 1048576^m * poch(1/2,m)^2 * prod(r=1,23, if(gcd(r,24)==1, poch(r/24,m)^2, 1)) / (poch(1,m)^6*poch(1/6,m)*poch(5/6,m)*poch(7/6,m)*poch(11/6,m)*poch(1/12,m)^2*poch(5/12,m)^2*poch(7/12,m)^2*poch(11/12,m)^2); }
read(concat(outdir,"c6_a2.seq"));  V6 = QV;
read(concat(outdir,"c12_a2.seq")); V12 = QV;
print("c6_a2  closed form matches for m<=60 : ", sum(m=0,60, if(qc6(m)==V6[m+1],0,1)) == 0);
print("c12_a2 closed form matches for m<=60 : ", sum(m=0,60, if(qc12(m)==V12[m+1],0,1)) == 0);
print("  qc6(3)  = ", qc6(3),  "   actual = ", V6[4]);
print("  qc12(3) = ", qc12(3), "   actual = ", V12[4]);
quit;
