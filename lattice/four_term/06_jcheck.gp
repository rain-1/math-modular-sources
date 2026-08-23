\r 06_jtest.gp
default(parisizemax, 8000000000);
/* control: Zagier A as a degenerate four-term row is not available; instead
   check the binomial-move image of Zagier A (rho=0, M=1, j=(0,1), a=4,c=1,d=-19,f=0,C=14)
   -- its local system is Zagier A's, so it must return deg J = 12. */
print("binom-move of Zagier A: ", jtest5(0,1,1,0,1, 4,1,-19,0,14));
quit;
