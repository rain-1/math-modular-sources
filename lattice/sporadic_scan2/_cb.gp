default(parisizemax, 4000000000);
read("identlib.gp");
default(realprecision, 60);
C = cuspbasis(12, 1, 50);
print("#C=", #C);
for(i=1,#C, if(abs(C[i][1] - 0.7372929961855962401764) < 1e-25, print("FOUND ", C[i][2], " ", C[i][1])));
quit
