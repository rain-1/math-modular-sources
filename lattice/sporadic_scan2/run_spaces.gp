{
for(i = 1, #LEVELS,
  my(N = LEVELS[i]);
  for(w = WLO, WHI,
    my(fn = Str("spaces/sp_", N, "_", w, ".txt"), t0, g, S, dim);
    if(#externstr(Str("test -f ", fn, " && echo y")) > 0, next);
    t0 = getabstime();
    g = ratgens(N, w);
    if(#g == 0,
      write(Str(fn,".part"), Str("N ", N, " w ", w, " DIM 0 PREC ", PREC));
    ,
      S = saturate(g);
      dim = matsize(S)[2];
      write(Str(fn,".part"), Str("N ", N, " w ", w, " DIM ", dim, " PREC ", PREC));
      for(c = 1, dim, write(Str(fn,".part"), Vec(S[,c])));
    );
    system(Str("mv ", fn, ".part ", fn));
    printf("N=%d w=%d  %dms\n", N, w, getabstime()-t0);
  );
);
}
quit
