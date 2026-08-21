x=1/7+O(5^25); test("planted 7x-1",[1,x],25,5);
y=3+2*5+4*5^3+O(5^25); test("random",[1,y],25,5); test("random quad",[1,y,y^2],25,5);
z=unitroot(-2,5,3,25); print("unitroot p=5: ",z," check ",z^2+2*z+125);
quit
