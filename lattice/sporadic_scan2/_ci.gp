default(parisizemax,4000000000);
read("identlib.gp");
default(realprecision,60);
C=cuspbasis(12,1,50);
print("n=",#C); for(i=1,min(6,#C), print(C[i][2]," ",C[i][1]));
C2=cuspbasis(8,1,50);
print("n8=",#C2); for(i=1,min(8,#C2), print(C2[i][2]," ",C2[i][1]));
quit
