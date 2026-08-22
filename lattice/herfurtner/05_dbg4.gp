read("05_jtest.gp");
chk(NP)={
 my(q,W,V,U,R,h=9,JSl,JJl,xq);
 JSl = jser(NP); JJl = x*JSl;
 q = nome(1,0,0,9,3,27,NP); xq = x/q;
 W = subst(JJl, x, q^h) * xq^h;
 V = x^9*(27*x^2-9*x+1);
 \\ solve for U = V*W / x^9 truncated to degree 12, then residual
 my(Wv = V*W/x^9);
 U = truncate(Wv + O(x^13));
 R = V*W - x^h*U;
 printf("NP=%d serprec(W)=%d ", NP, serprec(W,x));
 for(e=0,NP, if(polcoeff(R,e)!=0, printf("first nonzero at e=%d, val=%s\n",e,polcoeff(R,e)); return()));
 printf("all zero up to %d\n",NP);
};
chk(40); chk(56); chk(80);
quit;
