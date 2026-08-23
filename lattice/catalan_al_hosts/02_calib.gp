/* 02_calib.gp -- calibration: rebuild the level-8 (Zagier E) and level-16
   Catalan hosts from (F, Phi_0) inside M_1(N,chi_{-4}) x M_3(N,chi_{-4}). */
read("lattice/catalan_al_hosts/lib.gp");
NQ = 160;

{docase(nn, fv, pv) =
  my(mf1 = mfinit([nn,1,-4],4), mf3 = mfinit([nn,3,-4],4));
  my(fq = qexp(mf1, fv, NQ), pq = qexp(mf3, pv, NQ));
  print("  F  = ", vector(9,i,fq[i]));
  print("  P0 = ", vector(9,i,pq[i]));
  my(fs = vec2ser(fq), ps = vec2ser(pq));
  my(tq = Dinv(ps/fs));
  print("  t  = ", vector(9,i,polcoeff(tq,i-1)));
  my(qt = serreverse(tq));
  my(av = Vec(subst(fs,x,qt)));
  print("  A  = ", vector(9,i,av[i]));
  my(bv = Vec(subst(fs*Dinv2(ps), x, qt)));
  print("  B  = ", vector(7,i,bv[i]));
  av;
}
print("=== level 8: F = theta^2, Phi_0 = (1-8V_2)E ===");
mf1=mfinit([8,1,-4],4); print("  M1 basis:", vector(mfdim(mf1),i,mfcoefs(mfbasis(mf1)[i],6)));
mf3=mfinit([8,3,-4],4); print("  M3 basis:", vector(mfdim(mf3),i,mfcoefs(mfbasis(mf3)[i],6)));
quit;
