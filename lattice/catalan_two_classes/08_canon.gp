read("00_setup.gp");
/* canonical-source checks: Phi = F * D_q t  */
Dq(f) = my(r=0); for(n=0,NT-1, r += n*polcoeff(f,n)*q^n); r + O(q^NT);
c8  = F8*Dq(t8);   /* should be (1-8V2)E */
c16 = F16*Dq(x16); /* should be E        */
print("F8*Dq(t8)  - (1-8V2)E : ", c8  - ser(mkPhi(Sin,[[1,1],[2,-8]])) );
print("F16*Dq(x16) - E       : ", c16 - ser(mkPhi(Sin,[[1,1]])) );
print("t as a function of x (t = x(1+2x)/(1+4x)^2 ?): ", t8 - x16*(1+2*x16)/(1+4*x16)^2 );
