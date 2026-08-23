/* 01_jmap.gp -- recover the J-map U/V of the elliptic-K3 four-term row
 *   (n+1)^2 u_{n+1} = (11n^2+11n+4)u_n - (37n^2+3)u_{n-1} + 3(3n-1)(3n-2)u_{n-2}
 * class (rho;M,j1,j2) = (0;3,1,2), parameters (a,c,d,f,C) = (11,4,37,3,3).
 * Uses the machinery of ../four_term/06_jtest.gp and 13_jdetail.gp.
 */
\r /home/ubuntu/code/math-modular-sources/lattice/four_term/06_jtest.gp
\r /home/ubuntu/code/math-modular-sources/lattice/four_term/13_jdetail.gp
default(parisizemax, 8000000000);
v = jdetail(0,1,3,1,2, 11,4,37,3,3);
print();
ce = certify(0,1,3,1,2, 11,4,37,3,3, 4, -1, DMAX, NTERM);
uu = ce[4]; vv = ce[5];
print("UU := ", uu);
print("VV := ", vv);
print("content U = ", content(uu), "  content V = ", content(vv));
write("/home/ubuntu/code/math-modular-sources/lattice/k3_period/out/jmap.txt", "U = ", uu);
write("/home/ubuntu/code/math-modular-sources/lattice/k3_period/out/jmap.txt", "V = ", vv);
quit;
