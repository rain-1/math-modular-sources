/* 12_Gx.gp -- the q-product g(q) = G(x(q)) with G(x) = exp(int_0^x a(t) dt),
   a(x) = l(x)/(x sqrt(P) F).  Is G in Z[[x]]?  Is it algebraic / a power of
   a simple function of x?                                                       */
default(parisize, 4000000000);
read("lib.gp");
NM = 200;
{ for(k=1,3, my(R=ROWS[k], B=R[3], C=R[4], A=genrow(k,NM+2)[1], F, l, P, a, I, G);
  F = Ser(A,'x); l = 'x*Ser(vector(NM+3,n,A[n]/n),'x);
  P = 1 - 2*B*'x + (B^2-4*C)*'x^2 + O('x^(NM+3));
  a = l/('x*sqrt(P)*F);
  I = intformal(a);              /* int_0^x a dt */
  G = exp(I);
  print("=== ",NAM[k]);
  print("   a(x) = ",vector(6,j,polcoeff(a,j-1)));
  print("   G(x) = ",Ser(vector(11,j,polcoeff(G,j-1)),'x));
  my(bad=select(j->denominator(polcoeff(G,j))!=1, vector(NM,j,j)));
  print("   G in Z[[x]] for j<=",NM," ? ",if(#bad==0,"YES",Str("no, first bad j = ",bad[1..min(5,#bad)])));
  if(#bad, print("      denominators: ",vector(min(5,#bad),i,denominator(polcoeff(G,bad[i])))));
  /* algebraicity probe: is G a power of (1-lam1 x) or of P(x) ? */
  my(sq=sqrt(P), t1=log(G)/log(1-R[5]*'x+O('x^(NM+2))));
  print("   log G / log(1-lam1 x) = ",Ser(vector(6,j,polcoeff(t1,j-1)),'x));
  ); }
print();
print("=== cross-check: g(q) = G(x(q)) reproduces exp(sum c'(m)q^m/m) ? ===");
{ my(MM=40);
  for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4],A=genrow(k,MM+4)[1],F=Ser(A,'x),
       l='x*Ser(vector(MM+5,n,A[n]/n),'x), P=1-2*B*'x+(B^2-4*C)*'x^2+O('x^(MM+5)),
       G=exp(intformal(l/('x*sqrt(P)*F))), S=Setup(k,MM+4), xq=S[3], gq, cp, h);
    gq = subst(G,'x,xq);
    cp = vector(MM,m,polcoeff(S[4],m)/m);
    h = exp(Ser(concat([0],vector(MM,m,cp[m]/m)),'q));
    print("  ",NAM[k],": max difference of the first ",MM-4," coefficients = ",
      vecmax(vector(MM-4,n,abs(polcoeff(gq,n)-polcoeff(h,n)))))); }
quit;
