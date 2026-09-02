\\ 08_verify.gp -- PARI check of the two constant-term rules used in 07_supply.gp:
\\   (R1) a_0(Phi, infinity)      != 0  iff  psi_1 = 1
\\   (R2) a_0(Phi | [0,-1;1,0])   != 0  iff  psi_2 = 1
\\ and of the resulting dim(fold-regular) inside each nebentypus component.
\\ Entries live in different cyclotomic polmod rings, so they are embedded in C
\\ (consistently: t of Phi_n |-> exp(2 pi i/n)) before the rank is taken.
default(realprecision,60);
cycord(m) = {my(n,v=variable(m)); for(n=1,200, if(polcyclo(n,v)==m, return(n))); 0;}
toC(z) = {if(type(z)!="t_POLMOD", return(z*1.0)); my(m=z.mod, n=cycord(m)); if(n==0, error("not cyclotomic")); subst(lift(z), variable(m), exp(2*Pi*I/n));}
chk(N,k,me) =
{
  my(G=znstar(N,1), chi=znconreylog(G,me), mf, B, M, r, d, p);
  if(chareval(G,chi,-1)*2 != (k%2), print("N=",N," conrey=",me,": parity/weight mismatch, skipped"); return);
  mf = mfinit([N,k,[G,chi]],4);
  B = mfbasis(mf);
  d = #B;
  if(d==0, print("N=",N," k=",k," conrey=",me,"  dim 0"); return);
  M = matrix(2,d);
  for(j=1,d, M[1,j] = toC(mfcoefs(B[j],0)[1]));
  for(j=1,d, p=0; M[2,j] = toC(mfslashexpansion(mf,B[j],[0,-1;1,0],0,1,&p)[1]));
  r = matrank(M);
  print("N=",N," k=",k," conrey=",me," ord(eps)=",charorder(G,chi), "  dim(eps-comp)=",d,"  rank(a_0 at oo, a_0 at 0)=",r,"  dim(fold-regular)=",d-r);
}
chk(5,3,2);
chk(5,3,3);
chk(6,3,5);
chk(7,3,3);
chk(7,3,6);
chk(8,3,3);
chk(8,3,7);
chk(9,3,8);
chk(10,3,7);
chk(10,3,3);
chk(14,3,3);
chk(14,4,9);
chk(15,3,7);
chk(18,4,7);
quit;
