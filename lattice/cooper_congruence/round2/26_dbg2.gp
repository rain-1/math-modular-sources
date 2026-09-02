\\ 26_dbg2.gp -- independent check of all four rational forms of K mod p, to q^12000.
read("lib.gp");
G = [read("20c_gamma_s7.txt"), read("20c_gamma_s10.txt"), read("20c_gamma_s18.txt")];
N = #G[1];
print("N = ", N);
{
my(rows=[1,2,3,3], ps=[2,2,2,3],
   Ps=['q^11+'q^7+'q^5+'q,
       'q^15+'q^14+'q^13+'q^11+'q^5+'q^3+'q^2+'q,
       'q^7+'q^6+'q^2+'q,
       'q^3+'q^2+2*'q],
   Qs=['q^12+'q^10+'q^8+'q^6+'q^4+'q^2+1,
       'q^16+'q^12+'q^8+'q^4+1,
       'q^8+'q^4+1,
       2*'q^3+2]);
for(j=1,4,
  my(k=rows[j], p=ps[j], KK, S, bad);
  KK = sum(n=1,N, (G[k][n]%p)*'q^n) + O('q^(N+1));
  S  = Mod(1,p)*(Qs[j]*KK - Ps[j]);
  bad = sum(i=0,N, if(lift(polcoeff(S,i))!=0, 1, 0));
  print("  ",NAM[k]," p=",p,":  nonzero coefficients of Q*K-P in q^0..q^",N," : ", bad);
);
}
quit;
