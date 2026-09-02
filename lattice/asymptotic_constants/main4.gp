default(parisize,"12G");
read("lib.gp");
default(realprecision,400);
MQ = 1200;

Dop(S) = 'q*deriv(S,'q);
etaq(DL,EL)={
  my(r=1+O('q^(MQ+1)), s=0);
  for(i=1,#DL, s += DL[i]*EL[i]; for(m=1,MQ\DL[i], r *= (1-'q^(DL[i]*m)+O('q^(MQ+1)))^EL[i]));
  'q^(s/24)*r;
}

/* name, N, p2, p1, B, C, divisors, exponents */
ROWS = List();
listput(ROWS, ["alpha", 12, (2*n+1)*(10*n^2+10*n+4), -64*n^3, 10, 9, [1,2,3,4,6,12], [-4,4,4,-4,-4,4]]);
listput(ROWS, ["gamma",  6, (2*n+1)*(17*n^2+17*n+5), -1*n^3,  17, 72, [1,2,3,6], [-5,1,-1,5]]);
listput(ROWS, ["eps",    8, (2*n+1)*(12*n^2+12*n+4), -16*n^3, 12, 32, [1,2,4,8], [-4,2,-2,4]]);
listput(ROWS, ["zeta",   9, (2*n+1)*(9*n^2+9*n+3),   27*n^3,   9, 27, [1,3,9], [-3,0,3]]);
listput(ROWS, ["s7",     7, (2*n+1)*(13*n^2+13*n+4), 3*n*(9*n^2-1),    13, 49, [1,7], [-4,4]]);
listput(ROWS, ["s10",   10, 2*(2*n+1)*(3*n^2+3*n+1), 4*n*(16*n^2-1),    6, 25, [1,2,5,10], [-2,-2,2,2]]);
listput(ROWS, ["s18",   18, 2*(2*n+1)*(7*n^2+7*n+3), -12*n*(16*n^2-1), 14, 1, [1,2,3,6,9,18], [-6,6,12,-12,-6,6]]);

NN = 6000; MPTS = 20; HSTEP = 150;

{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], N=R[2], p2=R[3], p1=R[4], BB=R[5], CC=R[6], DL=R[7], EL=R[8]);
  print("################ ", nm, "  N=", N, "  B=", BB, "  C=", CC);
  my(uq = etaq(DL,EL));
  print("  u (eta quotient) = ", uq+O('q^8));
  my(xq = uq/(1+BB*uq+CC*uq^2), Fq = Dop(uq)/uq);
  my(qx = serreverse(xq+O('q^45)), Ax = subst(Fq+O('q^45),'q,qx));
  my(av = seqA((n+1)^3,p2,p1,44));
  print("  a_n from eta (x,F) vs recurrence, n<=40 : ", vector(41,j,polcoef(Ax,j-1)-av[j])==vector(41,j,0));
  print("  a_0..a_6 = ", vector(7,j,av[j]));
  /* exact recurrence to NN */
  my(A = seqA((n+1)^3,p2,p1,NN));
  print("  a_", NN, " integral: ", denominator(A[NN+1])==1, "   log10|a_N| = ", log(abs(A[NN+1]*1.))/log(10));
  my(sC=sqrt(CC*1.), l1=BB+2*sC, l2=BB-2*sC, xp=1/l1);
  my(Kn(m) = abs(A[m+1])*xp^m*m^(3/2));
  my(ns = vector(MPTS,j,NN-HSTEP*(j-1)));
  my(V = matrix(MPTS,MPTS,r,s,1/ns[r]^(s-1)), b = vector(MPTS,r,Kn(ns[r]))~);
  my(sol = matsolve(V,b), Kr = sol[1]);
  my(sol2 = matsolve(matrix(MPTS-2,MPTS-2,r,s,1/ns[r]^(s-1)), vector(MPTS-2,r,Kn(ns[r]))~));
  print("  K (Richardson ",MPTS," pts)  = ", Kr);
  print("  K (Richardson ",MPTS-2," pts)  = ", sol2[1]);
  print("  agreement digits          = ", -log(abs(Kr/sol2[1]-1))/log(10));
  my(Kcf = sqrt(N)/(2*Pi^(3/2))*sqrt(l1/(l1-l2)));
  print("  K (closed form)           = ", Kcf);
  print("  closed form vs Richardson : ", Kr/Kcf-1);
  print("  K*pi^{3/2}                = ", Kr*Pi^(3/2));
  print("  algdep(K pi^{3/2}, 4)     = ", algdep(Kr*Pi^(3/2),4));
  print("  algdep(K^2 pi^3, 2)       = ", algdep(Kr^2*Pi^3,2));
  /* fold formula, evaluated from the eta series at q_c */
  my(qc = exp(-2*Pi/sqrt(N)));
  my(uc = subst(truncate(uq),'q,qc), Du = subst(truncate(Dop(uq)),'q,qc), D2u = subst(truncate(Dop(Dop(uq))),'q,qc));
  print("  u(tau_c)                  = ", uc, "   1/sqrt(C) = ", 1/sC);
  print("  u(tau_c)-1/sqrt(C)        = ", uc-1/sC);
  my(Fv = Du/uc, DFv = D2u/uc - Fv^2);
  my(rr(z) = z/(1+BB*z+CC*z^2));
  my(r2 = -2*CC*uc/(1+BB*uc+CC*uc^2)^2);
  my(D2x = r2*Du^2, xpv = rr(uc));
  print("  x_+ = x(tau_c)            = ", xpv, "    1/lambda_1 = ", xp);
  print("  F(tau_c)                  = ", Fv);
  print("  DF(tau_c)                 = ", DFv);
  print("  DF - (sqrt(N)/(2 Pi)) F   = ", DFv - sqrt(N)/(2*Pi)*Fv);
  print("  (F/Du)(tau_c)             = ", Fv/Du, "    sqrt(C) = ", sC);
  print("  x''(h_c) = r''(u_c)       = ", r2, "   -2C^{3/2}/lam1^2 = ", -2*CC^(3/2)/l1^2);
  print("  D2x(tau_c)                = ", D2x);
  my(Kfold = sqrt(xpv*DFv^2/(2*Pi*(-D2x))));
  print("  K (fold formula)          = ", Kfold);
  print("  fold vs Richardson        = ", Kfold/Kr-1);
);
}
quit;
