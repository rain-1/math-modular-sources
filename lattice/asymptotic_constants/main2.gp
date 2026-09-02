default(parisize,"8G");
read("lib.gp");
NT = 300;

Dop(S) = 't*deriv(S,'t);

/* rows: [name, N, p2, p1, a, c]  with char poly lam^2-2a lam + c */
ROWS = List();
listput(ROWS, ["alpha", 12, (2*n+1)*(10*n^2+10*n+4), -64*n^3, 10, 64]);
listput(ROWS, ["gamma",  6, (2*n+1)*(17*n^2+17*n+5), -1*n^3,  17, 1]);
listput(ROWS, ["eps",    8, (2*n+1)*(12*n^2+12*n+4), -16*n^3, 12, 16]);
listput(ROWS, ["zeta",   9, (2*n+1)*(9*n^2+9*n+3),   27*n^3,   9, -27]);
listput(ROWS, ["s7",     7, (2*n+1)*(13*n^2+13*n+4), 3*n*(9*n^2-1),    13, -27]);
listput(ROWS, ["s10",   10, 2*(2*n+1)*(3*n^2+3*n+1), 4*n*(16*n^2-1),    6, -64]);
listput(ROWS, ["s18",   18, 2*(2*n+1)*(7*n^2+7*n+3), -12*n*(16*n^2-1), 14, 192]);

ratfit(G,U,d1,d2,NC)={
  my(cols=List(), Ms, ker);
  for(j=0,d1, listput(cols, Vec(truncate(U^j+O('t^(NC+1))),-(NC+1))));
  for(j=0,d2, listput(cols, Vec(truncate(-G*U^j+O('t^(NC+1))),-(NC+1))));
  Ms = Mat(vector(#cols,i,cols[i]~));
  ker = matker(Ms);
  ker;
}

{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], N=R[2], p2=R[3], p1=R[4], a=R[5], c=R[6]);
  my(M=mirror((n+1)^3,p2,p1,NT), tq=M[1], Fq=M[2]);
  my(CC=(a^2-c)/4, BB=a);
  my(P3 = 1-2*a*tq+c*tq^2);
  my(uq = 2*((1-a*tq)-sqrt(P3))/((a^2-c)*tq));
  my(G = Fq/Dop(uq));
  print("=== ", nm, "  N=",N,"  B=",BB,"  C=",CC);
  print("  u = ", uq+O('t^9));
  print("  u integral: ", denominator(content(truncate(uq)))==1);
  print("  check t - u/(1+Bu+Cu^2) = ", tq - uq/(1+BB*uq+CC*uq^2)+O('t^20));
  print("  G = F/Du = ", G+O('t^9));
  my(found=0);
  for(d=1,6,
    my(k1=ratfit(G,uq,d,d,4*d+30));
    if(#k1>0 && found==0, found=d; print("  rational fit degree ",d,": kernel dim ",#k1); print("    coeffs: ",k1));
  );
  if(found==0, print("  no rational fit up to degree 6"));
);
}
quit;
