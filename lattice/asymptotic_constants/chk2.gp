read("lib.gp");
print("alpha a_0..a_6 = ", vector(7,j,seqA((n+1)^3,(2*n+1)*(10*n^2+10*n+4),-64*n^3,8)[j]));
quit;
