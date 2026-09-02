default(parisize,"10G");
read("lib.gp");
NT = 300;
Dop(S) = 't*deriv(S,'t);

ROWS = List();
listput(ROWS, ["alpha", 12, (2*n+1)*(10*n^2+10*n+4), -64*n^3, 10, 64]);
listput(ROWS, ["gamma",  6, (2*n+1)*(17*n^2+17*n+5), -1*n^3,  17, 1]);
listput(ROWS, ["eps",    8, (2*n+1)*(12*n^2+12*n+4), -16*n^3, 12, 16]);
listput(ROWS, ["zeta",   9, (2*n+1)*(9*n^2+9*n+3),   27*n^3,   9, -27]);
listput(ROWS, ["s7",     7, (2*n+1)*(13*n^2+13*n+4), 3*n*(9*n^2-1),    13, -27]);
listput(ROWS, ["s10",   10, 2*(2*n+1)*(3*n^2+3*n+1), 4*n*(16*n^2-1),    6, -64]);
listput(ROWS, ["s18",   18, 2*(2*n+1)*(7*n^2+7*n+3), -12*n*(16*n^2-1), 14, 192]);

/* log of eta-quotient factor prod_n (1-q^{dn}) as a t-series to order NC */
lgpr(d,NC)={
  my(s=O('t^(NC+1)));
  for(nn=1,NC\d, s += log(1-'t^(d*nn)+O('t^(NC+1))));
  s;
}

etaid(uq,DIVS,NC)={
  my(L=Vec(truncate(log(uq/'t+O('t^(NC+1)))),-(NC+1)), Ms, rhs, sol);
  L = vector(NC,k,L[k+1]);
  Ms = Mat(vector(#DIVS,i,my(v=Vec(truncate(lgpr(DIVS[i],NC)),-(NC+1))); vector(NC,k,v[k+1])~));
  sol = matsolve(Ms[1..#DIVS,]~*0,[]);
  0;
}

{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], N=R[2], p2=R[3], p1=R[4], a=R[5], c=R[6]);
  my(M=mirror((n+1)^3,p2,p1,NT), tq=M[1], Fq=M[2]);
  my(CC=(a^2-c)/4, BB=a);
  my(uq = 2*((1-a*tq)-sqrt(1-2*a*tq+c*tq^2))/((a^2-c)*tq));
  print("=== ", nm, "  N=",N,"  B=",BB,"  C=",CC);
  print("  u integral                : ", denominator(content(truncate(uq)))==1);
  print("  t - u/(1+Bu+Cu^2)         : ", tq - uq/(1+BB*uq+CC*uq^2));
  print("  F - Du/u                  : ", Fq - Dop(uq)/uq);
  print("  P3(t) - ((1-Cu^2)/(1+Bu+Cu^2))^2 : ", (1-2*a*tq+c*tq^2) - ((1-CC*uq^2)/(1+BB*uq+CC*uq^2))^2);
  print("  u = ", uq+O('t^13));
  write("useries.txt", [nm, N, BB, CC, Vec(truncate(uq),-201)]);
);
}
quit;
