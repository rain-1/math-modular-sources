read("05_jtest.gp");
chk(gam,h)={
 my(q,W,V,U,R,NP=56,JSl,JJl,xq);
 JSl = jser(NP); JJl = x*JSl;
 q = nome(1,0,0,9,3,27,NP); xq = x/q;
 W = subst(JJl, x, gam*q^h) * xq^h / gam;
 V = x^9*(27*x^2-9*x+1);
 U = truncate(V*W/x^9 + O(x^13));
 R = V*W - x^h*U;
 printf("gam=%s h=%d : ", gam, h);
 for(e=0,50, if(polcoeff(R,e)!=0, printf("first nonzero at e=%d\n",e); return()));
 printf("OK, U = %s\n",U);
};
chk(1,9); chk(-1,9);
quit;
