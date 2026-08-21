{
for(i = 1, #LEVELS,
  my(N = LEVELS[i]);
  for(w = 1, 3,
    my(CL = charorbits(N, w));
    for(ci = 1, #CL,
      my(m = CL[ci][1]);
      my(fn = Str("spaces/sp_", N, "_", w, "_", m, ".txt"), out = Str("prep/pr_", N, "_", w, "_", m, ".txt"));
      if(#externstr(Str("test -f ", fn, " && echo y")) == 0, next);
      if(#externstr(Str("test -f ", out, " && echo y")) > 0, next);
      iferr(doprep(fn, Str(out, ".part")), E, print("ERR ", N, " ", w, " ", m, " ", E));
      system(Str("mv ", out, ".part ", out, " 2>/dev/null"));
    );
  );
);
}
quit
