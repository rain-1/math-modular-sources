read("e2.gp");
default(realprecision,60);
\\ sanity: E2*(i) should be 3/pi ; E2*(rho)=0
print("E2*(i)    = ", E2star(I), "   3/Pi = ", 3/Pi);
print("E2*(rho)  = ", E2star((1+I*sqrt(3))/2));
print("E2*(0.03+0.011I) vs modular check:");
t = 0.03+0.011*I;
print(" direct  ", E2star(t));
print(" via -1/t ", E2star(-1/t)/t^2);
quit;
