default(realprecision, 70);
xi = -0.00050455459344136708862542545797071516117215446057914345028836870400436268;
print("xi = ", xi);
{cands = [["1",1],["zeta(3)",zeta(3)],["pi^2",Pi^2],["pi^4",Pi^4],["zeta(4)",zeta(4)],
         ["G",Catalan],["zeta(3)*Pi",zeta(3)*Pi],["zeta(2)^2",zeta(2)^2],
         ["L(-3,2)",lfun(-3,2)],["L(-3,3)",lfun(-3,3)],["L(-4,3)",lfun(-4,3)],
         ["L(-4,4)",lfun(-4,4)],["L(-3,4)",lfun(-3,4)],["L(-8,4)",lfun(-8,4)],
         ["L(-7,4)",lfun(-7,4)],["L(13,4)",lfun(13,4)],["L(17,4)",lfun(17,4)],
         ["L(-52,4)",lfun(-52,4)],["L(-68,4)",lfun(-68,4)],["L(-104,4)",lfun(-104,4)],
         ["L(-4,2)",lfun(-4,2)],["L(13,2)",lfun(13,2)],["L(17,2)",lfun(17,2)],
         ["L(13,3)",lfun(13,3)],["L(17,3)",lfun(17,3)],["L(-51,4)",lfun(-51,4)],
         ["L(-119,4)",lfun(-119,4)],["L(-8,3)",lfun(-8,3)],["L(8,4)",lfun(8,4)],
         ["L(5,4)",lfun(5,4)],["L(-3,5)",lfun(-3,5)],["zeta(5)",zeta(5)]];}
{for(j=1,#cands,
  my(r = xi/cands[j][2], q = bestappr(r, 10^10));
  if(q!=0 && abs(r-q) < abs(r)*1e-50,
     print("  HIT: xi = ", q, " * ", cands[j][1])));}
print("\n-- two-term lindep against pairs with small height --");
{for(j=1,#cands, for(l=j+1,#cands,
  my(v = lindep([xi, cands[j][2], cands[l][2]], 55));
  if(v!=0 && vecmax(abs(v))<10^7 && v[1]!=0,
     print("  ", cands[j][1], " & ", cands[l][1], " : ", v~))));}
quit
