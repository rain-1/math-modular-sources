read("05_jtest.gp");
dbg(M,j1,j2,A,B,C,dmax,N)={
 my(q,Jh,Jt,W,mat,ker,neq);
 q = nome(M,j1,j2,A,B,C,N);
 for(h=1,12,
   Jh = subst(JS, x, x^h);
   Jt = subst(Jh, x, q);
   W  = Jt*x^h;
   neq = 2*dmax+6; if(neq > serprec(W,x)-4, neq = serprec(W,x)-4);
   mat = matrix(neq, 2*dmax+2);
   for(i=1, neq, my(e=i-1);
     for(k=0,dmax, mat[i,k+1]=if(e-k>=0, polcoeff(W,e-k),0));
     for(k=0,dmax, mat[i,dmax+2+k]=if(e==h+k,-1,0)));
   ker = matker(mat);
   printf("h=%d serprec=%d neq=%d dimker=%d\n",h,serprec(W,x),neq,#ker);
 );
};
print("--- Zagier B (9,3,27) ---"); dbg(1,0,0,9,3,27,12,NS);
quit;
