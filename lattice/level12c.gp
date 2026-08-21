read("level12b.gp");
print("ratios: ",vector(6,i,A[M-i]*1./A[M-i-1]));
{for(ord=2,3,for(d=2,8, my(rows=[]); for(n=ord,M-1, my(r=[]); for(j=0,ord, for(e=0,d, r=concat(r,[n^e*A[n+1-j]]))); rows=concat(rows,[r]));
 my(Mx=matrix(#rows,(ord+1)*(d+1),i,j,rows[i][j])); print("ord ",ord," deg ",d," kerdim ",#matker(Mx))));}
\q
