default(parisize,"8G");
read("lib.gp");
default(realprecision,200);
NT = 600;

ev(S,z) = subst(truncate(S),'t,z);
Dop(S) = 't*deriv(S,'t);

/* rows: [name, N, p2, p1, Ppoly(x) whose smallest root is x_+ ] */
ROWS = List();
listput(ROWS, ["alpha", 12, (2*n+1)*(10*n^2+10*n+4), -64*n^3, 1-20*x+64*x^2]);
listput(ROWS, ["gamma",  6, (2*n+1)*(17*n^2+17*n+5), -1*n^3,  1-34*x+1*x^2]);
listput(ROWS, ["eps",    8, (2*n+1)*(12*n^2+12*n+4), -16*n^3, 1-24*x+16*x^2]);
listput(ROWS, ["zeta",   9, (2*n+1)*(9*n^2+9*n+3),   27*n^3,  1-18*x-27*x^2]);
listput(ROWS, ["s7",     7, (2*n+1)*(13*n^2+13*n+4), 3*n*(9*n^2-1),    1-26*x-27*x^2]);
listput(ROWS, ["s10",   10, 2*(2*n+1)*(3*n^2+3*n+1), 4*n*(16*n^2-1),   1-12*x-64*x^2]);
listput(ROWS, ["s18",   18, 2*(2*n+1)*(7*n^2+7*n+3), -12*n*(16*n^2-1), 1-28*x+192*x^2]);

{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], N=R[2], p2=R[3], p1=R[4], PP=R[5]);
  my(M=mirror((n+1)^3,p2,p1,NT), tq=M[1], Fq=M[2]);
  my(tt=truncate(tq), FF=truncate(Fq));
  print("=== ", nm, "   N=", N);
  print("  t integral: ", denominator(content(tt))==1, "   F integral: ", denominator(content(FF))==1);
  print("  t = ", tq+O('t^7));
  print("  F = ", Fq+O('t^7));
  my(qc0 = exp(-2*Pi/sqrt(N)));
  my(dt=deriv(tt), ddt=deriv(dt), qc=qc0);
  for(it=1,60, qc = qc - subst(dt,'t,qc)/subst(ddt,'t,qc));
  print("  q_c (Newton on t'=0)      = ", qc);
  print("  exp(-2 Pi/sqrt(N))        = ", qc0);
  print("  q_c/exp(-2Pi/sqrtN) - 1   = ", qc/qc0-1);
  my(xp = vecsort(polroots(PP),abs)[1]);
  print("  x_+ (smallest root of P)  = ", xp);
  print("  t(q_c) - x_+              = ", ev(tq,qc0)-xp);
  print("  tail size |c_NT q_c^NT|   = ", abs(polcoef(tt,NT))*qc0^NT);
  my(D1t=ev(Dop(tq),qc0), D2t=ev(Dop(Dop(tq)),qc0));
  my(Fv=ev(Fq,qc0), DF=ev(Dop(Fq),qc0));
  print("  Dt(tau_c)                 = ", D1t);
  print("  D2t(tau_c)                = ", D2t);
  print("  F(tau_c)                  = ", Fv);
  print("  DF(tau_c)                 = ", DF);
  print("  DF/F                      = ", DF/Fv, "   (k sqrt(N)/(4 Pi), k=2: ", 2*sqrt(N)/(4*Pi), ")");
  my(K2 = xp*DF^2/(2*Pi*(-D2t)));
  print("  K (fold formula)          = ", sqrt(K2));
  write("fold_out.txt", [nm, N, qc0, xp, Fv, DF, D2t, sqrt(K2)]);
);
}
quit;
