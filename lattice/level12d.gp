read("level12b.gp");
{for(ord=4,10,for(d=2,3, my(rows=[]); for(n=ord,M-1, my(r=[]); for(j=0,ord, for(e=0,d, r=concat(r,[n^e*A[n+1-j]]))); rows=concat(rows,[r]));
 my(Mx=matrix(#rows,(ord+1)*(d+1),i,j,rows[i][j]),K=matker(Mx)); print("ord ",ord," deg ",d," kerdim ",#K);
 if(#K>0, my(v=K[,1]); v=v/content(v); print(vector(ord+1,j,Polrev(vector(d+1,e,v[(j-1)*(d+1)+e]),n))); break(2))));}
\q
