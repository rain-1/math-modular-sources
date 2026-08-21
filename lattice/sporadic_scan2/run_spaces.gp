{
for(i = 1, #LEVELS,
  my(N = LEVELS[i], psi = if(N==1, 1, N*prod(k=1,#factor(N)~, 1 + 1/factor(N)[k,1])));
  for(w = WLO, WHI,
    my(CL = charorbits(N, w), sturm = ceil(w*psi/12));
    for(ci = 1, #CL,
      my(m = CL[ci][1], o = CL[ci][2], ph = eulerphi(o));
      if(ph*(sturm+2) + 10 > PREC, next);
      my(fn = Str(SPDIR, "/sp_", N, "_", w, "_", m, ".txt"), t0, g, S, dim);
      if(#externstr(Str("test -f ", fn, " && echo y")) > 0, next);
      t0 = getabstime();
      g = orbgens(N, w, m);
      if(#g == 0,
        write(Str(fn,".part"), Str("N ", N, " w ", w, " CHI ", m, " ORD ", o, " DIM 0 PREC ", PREC));
      ,
        S = if(PREC > 260, saturate_fast(g, 210), saturate(g));
        dim = matsize(S)[2];
        write(Str(fn,".part"), Str("N ", N, " w ", w, " CHI ", m, " ORD ", o, " DIM ", dim, " PREC ", PREC));
        for(c = 1, dim, write(Str(fn,".part"), Vec(S[,c])));
      );
      system(Str("mv ", fn, ".part ", fn));
      printf("N=%d w=%d chi=%d ord=%d %dms\n", N, w, m, o, getabstime()-t0);
    );
  );
);
}
quit
