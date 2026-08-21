default(parisizemax, 4000000000);
read("identlib.gp");
default(realprecision, 60);
C = cuspbasis(12, 1, 50);
print("#C=", #C);
for(i=1,#C, if(abs(C[i][1] - 0.73729299618) < 1e-8, print("NEAR ", C[i][2], " ", C[i][1])));
for(i=1,20, print(C[i][2], " = ", C[i][1]));
quit
