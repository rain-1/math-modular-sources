read("lib.gp");
M = 300;
{
for(k=1,3,
  my(S=Setup(k,M+4), u=S[1], Ph=S[4], c, flat, PF, r, mat, ker, DD, v, P, Q, res);
  c = vector(M,m,polcoeff(Ph,m));
  flat = vector(M);
  for(n=1,M, my(s=0); fordiv(n,e, s += moebius(e)*psin(k,e)*e*c[n/e]); flat[n]=s);
  PF = sum(n=1,M, flat[n]*'q^n) + O('q^(M+1));
  r = PF/(Ph+O('q^(M+1)));
  DD = 6;
  mat = matrix(60, 2*DD+2, i, j,
     if(j<=DD+1, polcoeff(u^(j-1)+O('q^(M+1)), i-1), -polcoeff(r*u^(j-DD-2)+O('q^(M+1)), i-1)));
  ker = matker(mat);
  print("row ",NAM[k],"  dim ker (60 eqns, ",2*DD+2," unknowns) = ",#ker);
  if(#ker>0,
    v = ker[,1];
    P = sum(i=0,DD, v[i+1]*u^i);
    Q = sum(j=0,DD, v[DD+2+j]*u^j);
    res = P - r*Q;
    print("   coeffs of P: ", vector(DD+1,i,v[i]), "   coeffs of Q: ", vector(DD+1,i,v[DD+1+i]));
    print("   max q-order checked: 250;  first nonzero coeff of P - r*Q at order: ",
      valuation(res+O('q^250),'q))
  );
);
}
quit;
