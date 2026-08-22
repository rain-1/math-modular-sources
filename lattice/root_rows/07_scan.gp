read("lattice/root_rows/lib.gp");
p = 2^61-1;
readrow(name)={readvec(Str("lattice/root_rows/rows_",name,".txt"))};

\\ mod-p kernel dimension of the (order,degree) ansatz, with margin equations
scan(a, maxo, maxd)={
 my(N=#a-1);
 print("   ord\\deg  ", vector(maxd+1,j,j-1));
 for(o=1,maxo,
   my(line=vector(maxd+1));
   for(d=0,maxd,
     my(nu=(o+1)*(d+1));
     if(nu+20 > N-o, line[d+1]=-1,
       my(rows=[]);
       for(n=o+1,N, my(r=vector(nu)); my(c=1);
         for(j=0,o, for(e=0,d, r[c]=Mod(n,p)^e*Mod(a[n+1-j],p); c++));
         rows=concat(rows,[r]));
       my(Mx=matrix(#rows,nu,i,j,rows[i][j]));
       line[d+1]=nu-matrank(Mx)));
   print("   ",o,"      ",line));
};
