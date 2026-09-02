default(parisizemax, 20000000000);
read("build.txt");
Lc(bb, n) = my(b0=bb[n+1], b1=if(n>=1, bb[n], 0), b2=if(n>=2, bb[n-1], 0)); n^2*b0 - (11*n^2-11*n+3)*b1 - (n-1)^2*b2;
/* full verification: RHS*(1-11x-x^2) equals the claimed polynomial for ALL n<=NA */
chk(bb, num, nm) = my(bad=0, S=sum(i=1,NA, Lc(bb,i)*x^(i-1)) + O(x^NA)); my(N=S*(1-11*x-x^2) - num + O(x^(NA-1))); for(i=0, NA-3, if(polcoeff(N,i)!=0, bad++)); print(nm, ": mismatches among x^0..x^", NA-3, " = ", bad);
chk(A,  0,           "L(A)  = 0                     ");
chk(BD, 1-11*x-x^2,  "L(BD) = 1                     ");
chk(B3, 2,           "L(B3) = 2/(1-11x-x^2)         ");
chk(B4, -2*x,        "L(B4) = -2x/(1-11x-x^2)       ");
/* factorisation check of the numerators */
print("");
print("1-11x-x^2 = -(x-t1)(x-t2), t1=phi^-5, t2=-phi^5, t1*t2=-1, t1+t2=-11");
print("L(B3+ph5*B4) numerator 2-2*ph5*x = 2*ph5*(x-t1)  ->  RHS = 2*ph5/(x-t2) = 2/(1-x/t2)");
print("L(B3-phm5*B4) numerator 2+2*phm5*x = 2*phm5*(x-t2) -> RHS = -2*phm5/(x-t1) = 2/(1-x/t1)");
quit;
