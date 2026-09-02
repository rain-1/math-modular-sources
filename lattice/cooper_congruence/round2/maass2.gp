\\ maass2.gp -- robust fhat: closed form (E_2^*, E_4 by SL_2(Z) reduction) everywhere,
\\ q-series fallback exactly at the zeros of F (the polar CM points).
read("maass.gp");
FSER = vector(3);
{ initfser(k,NQ) = FSER[k] = mkfser(k,NQ); }
{ fhatR(k,t) = my(N=LEV[k], Fv0, A);
  Fv0 = Fmod(k,t);
  if(abs(Fv0) < 1e-25, return(fhatQ(t,N,FSER[k],#Vec(FSER[k][1])-3)));
  fhatC(k,t)[2];
}
