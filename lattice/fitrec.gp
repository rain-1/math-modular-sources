fitrec(A,ord,d)={my(M=#A,rows=[]);
 for(n=ord,M-1, my(r=[]); for(j=0,ord, for(e=0,d, r=concat(r,[n^e*A[n+1-j]]))); rows=concat(rows,[r]));
 my(Mx=matrix(#rows,(ord+1)*(d+1),i,j,rows[i][j]),K=matker(Mx));
 if(#K==0,return(0)); my(v=K[,1]); v=v/content(v);
 vector(ord+1,j,Polrev(vector(d+1,e,v[(j-1)*(d+1)+e]),n))};
