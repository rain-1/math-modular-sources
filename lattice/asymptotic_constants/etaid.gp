default(parisize,"8G");
read("lib.gp");
NT = 120;
NC = 60;
Dop(S) = 't*deriv(S,'t);

ROWS = List();
listput(ROWS, ["alpha", 12, (2*n+1)*(10*n^2+10*n+4), -64*n^3, 10, 64]);
listput(ROWS, ["gamma",  6, (2*n+1)*(17*n^2+17*n+5), -1*n^3,  17, 1]);
listput(ROWS, ["eps",    8, (2*n+1)*(12*n^2+12*n+4), -16*n^3, 12, 16]);
listput(ROWS, ["zeta",   9, (2*n+1)*(9*n^2+9*n+3),   27*n^3,   9, -27]);
listput(ROWS, ["s7",     7, (2*n+1)*(13*n^2+13*n+4), 3*n*(9*n^2-1),    13, -27]);
listput(ROWS, ["s10",   10, 2*(2*n+1)*(3*n^2+3*n+1), 4*n*(16*n^2-1),    6, -64]);
listput(ROWS, ["s18",   18, 2*(2*n+1)*(7*n^2+7*n+3), -12*n*(16*n^2-1), 14, 192]);

lgpr(d)={
  my(s=O('t^(NC+1)));
  for(m=1,NC\d, s += log(1-'t^(d*m)+O('t^(NC+1))));
  s;
}

tryset(uq,DIVS)={
  my(L=Vec(truncate(log(uq/'t+O('t^(NC+1)))),-(NC+1)), Ms, sol);
  L = vector(NC,k,L[k+1])~;
  Ms = Mat(vector(#DIVS,i,my(v=Vec(truncate(lgpr(DIVS[i])),-(NC+1))); vector(NC,k,v[k+1])~));
  sol = matinverseimage(Ms,L);
  sol;
}

{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], N=R[2], p2=R[3], p1=R[4], a=R[5], c=R[6]);
  my(M=mirror((n+1)^3,p2,p1,NT), tq=M[1]);
  my(uq = 2*((1-a*tq)-sqrt(1-2*a*tq+c*tq^2))/((a^2-c)*tq));
  print("=== ", nm, " N=", N);
  foreach([N,2*N,3*N,4*N], LV,
    my(DIVS=divisors(LV), s=tryset(uq,DIVS));
    if(type(s)=="t_COL" && #s>0,
      print("  eta quotient on divisors of ",LV,": ",DIVS," exponents ",s~,"  sum=",vecsum(Vec(s)),"  ord_inf=",sum(j=1,#DIVS,DIVS[j]*s[j])/24);
      break;
    );
  );
);
}
quit;
