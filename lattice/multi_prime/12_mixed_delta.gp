\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/12_mixed_delta.gp
\\ Master formula for the two-slope-prime cusp-move rows D.2, D.3, D.4, D.6
\\ (sqrt(Domb) orbit).  Measured: sigma_2 = 4 (D.2, D.6), 2 (D.3, D.4); sigma_3 = 1 for all.
\\ k = 2 for all four (CUSP_MOVE_PROGRAM section 4.1), Lambda = 16, |lambda_2| = 12 resp. 4.
\\ None of them decays, so F > 0 is automatic and no irrationality content exists;
\\ the number below measures only what the second prime is worth.
default(parisizemax, 8000000000);
default(realprecision, 38);

\\ [name, k, Lambda, |lambda_2|, sigma_2, sigma_3]
{RW = [["D.2", 2, 16, 12, 4, 1], ["D.3", 2, 16, 4, 2, 1],
       ["D.4", 2, 16, 4, 2, 1], ["D.6", 2, 16, 12, 4, 1]];}

dsc(id, ie) = {
  my(dd = RW[id], ee = RW[ie], best = vector(4), br = vector(4));
  for(j = 1, 400,
    my(rr = j*0.01, kk = 2*max(1, rr), g2, g3, gs, ff, hh, dl);
    g2 = min(ee[5]*rr, dd[5])*log(2);
    g3 = min(ee[6]*rr, dd[6])*log(3);
    gs = [0.0, g3, g2, g2+g3];
    for(w = 1, 4,
      ff = 0.5*(kk + rr*log(ee[4]) + log(dd[4]) - gs[w]);
      hh = ff + log(dd[3]) - log(dd[4]);
      dl = if(hh <= 0, -1.0, 1 - ff/hh);
      if(dl > best[w], best[w] = dl; br[w] = rr)));
  printf("  decayer %-4s engine %-4s :  none %6.4f | p=3 %6.4f | p=2 %6.4f | BOTH %6.4f  (r=%4.2f)\n",
    dd[1], ee[1], best[1], best[2], best[3], best[4], br[4]);
  best[4];
};

print("### master formula on the two-slope-prime sqrt(Domb) rows");
print("    (all four have |lambda_2| > 1: the linear forms GROW, F > 0 always,");
print("     so these numbers are worthiness bookkeeping, not irrationality)");
{for(i = 1, 4, for(j = 1, 4, if(i != j, dsc(i, j))));}
quit
