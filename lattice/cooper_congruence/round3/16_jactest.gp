\\ 16_jactest.gp -- sanity checks on the weak Jacobi generators.
default(parisize, 4000000000);
NQTEST = 1;
read("15_jaclib.gp");
P21 = mkphi21();
P01 = mkphi01();
P12 = mkphi12();
print("phi_{-2,1}:  q^0 coefficients r=-3..3: ", vector(7,i,jcoef(P21,0,i-4)));
print("             q^1 coefficients r=-3..3: ", vector(7,i,jcoef(P21,1,i-4)));
print("phi_{0,1}:   q^0 coefficients r=-3..3: ", vector(7,i,jcoef(P01,0,i-4)));
print("             q^1 coefficients r=-3..3: ", vector(7,i,jcoef(P01,1,i-4)));
print("phi_{-1,2}:  q^0 coefficients r=-3..3: ", vector(7,i,jcoef(P12,0,i-4)));
print("             q^1 coefficients r=-4..4: ", vector(9,i,jcoef(P12,1,i-5)));
print("E4: ", vector(5,i,polcoeff(E4S,2*(i-1))));
print("Delta: ", vector(5,i,polcoeff(DELS,2*(i-1))));
\\ a weak Jacobi form of weight k index m has c(n,r) depending only on 4mn-r^2 and r mod 2m
{ chkindex(A,m,name) = my(bad=0);
   for(n=0,6, for(r=-2*m-4, 2*m+4,
     my(D=4*m*n-r^2, c=jcoef(A,n,r));
     for(n2=0,6, for(r2=-2*m-4,2*m+4,
       if(4*m*n2-r2^2==D && (r-r2)%(2*m)==0 && jcoef(A,n2,r2)!=c, bad++)))));
   print(name, ": index-", m, " consistency violations = ", bad); }
chkindex(P21,1,"phi_{-2,1}");
chkindex(P01,1,"phi_{0,1}");
chkindex(P12,2,"phi_{-1,2}");
quit;
