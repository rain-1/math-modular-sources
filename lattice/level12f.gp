read("level12b.gp");
ord=7;d=3; rows=[]; for(n=ord,M-1, r=[]; for(j=0,ord, for(e=0,d, r=concat(r,[n^e*A[n+1-j]]))); rows=concat(rows,[r]));
Mx=matrix(#rows,(ord+1)*(d+1),i,j,rows[i][j]); K=matker(Mx);
print("kernel dim ",#K);
\\ try to find a degree-2 combination in the kernel: look for v in span with n^3 coefficients zero
P=vector(#K,k,vector(ord+1,j,Polrev(vector(d+1,e,K[(j-1)*(d+1)+e,k]),n)));
\\ eliminate n^3 between the two kernel vectors
c1=polcoeff(P[1][1],3); c2=polcoeff(P[2][1],3);
W=vector(ord+1,j,c2*P[1][j]-c1*P[2][j]); W=W/content(W);
print("combo: ",W);
print("max degree: ",vecmax(vector(ord+1,j,poldegree(W[j]))));
\q
