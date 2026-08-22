read("05_jtest.gp");
{
 my(q,Jh,Jt,W,mat,ker,neq,dmax=12,h=9,V,U,R);
 q = nome(1,0,0,9,3,27,NS);
 Jh = subst(JS, x, x^h); Jt = subst(Jh, x, q); W = Jt*x^h;
 printf("serprec(W)=%d val(W)=%d\n", serprec(W,x), valuation(W,x));
 neq = 30;
 mat = matrix(neq, 2*dmax+2);
 for(i=1,neq, my(e=i-1);
   for(k=0,dmax, mat[i,k+1]=if(e-k>=0, polcoeff(W,e-k),0));
   for(k=0,dmax, mat[i,dmax+2+k]=if(e==h+k,-1,0)));
 ker = matker(mat);
 printf("dimker=%d\n",#ker);
 V = sum(k=0,dmax, ker[,1][k+1]*x^k);
 U = sum(k=0,dmax, ker[,1][dmax+2+k]*x^k);
 printf("V = %s\n", V); printf("U = %s\n", U);
 R = V*W - x^h*U;
 printf("serprec(R)=%d\n", serprec(R,x));
 for(e=0,45, if(polcoeff(R,e)!=0, printf("first nonzero coeff at e=%d\n",e); break));
}
quit;
