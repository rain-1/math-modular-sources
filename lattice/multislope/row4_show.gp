default(parisizemax, 6000000000);
ciB = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recB.txt");
QB = vector(36, k, my(i=k-1); sum(j=0,7, ciB[i*8+j+1]*'n^j));
for(k=1,4, print("Q_",k-1," = ", QB[k]));
print("Q_34 = ", QB[35]);
print("lead coeffs (coeff of n^7 in Q_i), i=0..35: ", vector(36,k,polcoeff(QB[k],7,'n)));
quit
