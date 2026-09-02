/* first coefficients of BD and Bnew = B3+phi^5 B4 = (u+v*sqrt5)/2, u=2B3+11B4, v=5B4 */
read("build.txt");
print("n : BD_n");
for(n=0,9, print("  ",n," : ", BD[n+1]));
print("");
print("n : B_new,n = a + b*sqrt5   (a=B3+11/2*B4, b=5/2*B4)");
for(n=0,9, print("  ",n," : a = ", B3[n+1]+11/2*B4[n+1], " ,  b = ", 5/2*B4[n+1]));
print("");
print("n : B_new',n = a' + b'*sqrt5 = B3 - phi^-5 B4  (a'=B3+11/2*B4, b'=-5/2*B4)");
for(n=0,9, print("  ",n," : a' = ", B3[n+1]+11/2*B4[n+1], " ,  b' = ", -5/2*B4[n+1]));
print("");
print("B3 = ", vector(10,i,B3[i]));
print("B4 = ", vector(10,i,B4[i]));
quit;
