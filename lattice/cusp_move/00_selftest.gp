/* 00_selftest.gp -- sanity checks of the cusp-move machinery on Zagier A,C,F. */
read("lib.gp");
default(parisize, 2^29);

zagA = [7*nv^2+7*nv+2, -8*nv^2, "Zagier A"];
zagC = [10*nv^2+10*nv+3, 9*nv^2, "Zagier C"];
zagF = [17*nv^2+17*nv+6, 72*nv^2, "Zagier F"];
N = 40;

print("A: ", Vec(seqA(zagA,6)), "  C: ", Vec(seqA(zagC,6)), "  F: ", Vec(seqA(zagF,6)));
aC = seqA(zagC, N);

{
for(i = 1, 2,
  my(lam, w, rw);
  lam = if(i==1, 1, 9);
  w = cmove(aC, lam, 1);
  rw = fitrow(w, 2);
  if(rw == 0, print("  fit FAILED for lam=", lam),
     print("  lam=", lam, " -> P=", rw[1], "  Q=", rw[2],
           "  verify=", checkrow(rw, w),
           "  (a,b,d)=(", rowa(rw), ",", rowb(rw), ",", rowd(rw), ")"));
);
}
print("  C -(1)-> ", Vec(cmove(aC,1,1)[1..7]), "   Zagier A: ", Vec(seqA(zagA,6)));
print("  C -(9)-> ", Vec(cmove(aC,9,1)[1..7]), "   (-1)^n F: ", Vec(vector(7,i,(-1)^(i-1)*seqA(zagF,6)[i])));

bC = seqB(zagC, N);
{
for(i = 1, 2,
  my(lam, bh, bt, w, rw, btrue, agree, bn, ag2);
  lam = if(i==1, 1, 9);
  bh = seqBg(zagC, N, (k) -> lam^(k-1));
  bt = cmove(bh, lam, 1);
  w  = cmove(aC, lam, 1);
  rw = fitrow(w, 2);
  btrue = seqB(rw, N);
  agree = 1;
  for(n = 0, N, if(bt[n+1] != btrue[n+1], agree = 0; break));
  print("  lam=", lam, " : B^# = cmove(Bhat) -> ", if(agree, Str("EXACT to n=",N), "FAILS"));
  bn = cmove(bC, lam, 1); ag2 = 1;
  for(n = 0, N, if(bn[n+1] != btrue[n+1], ag2 = 0; break));
  print("            naive B^# = cmove(B)   -> ", if(ag2, "EXACT", "FAILS (expected)"));
);
}
quit;
