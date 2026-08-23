/* (g) recurrence fit for the central pair, over Q, via matker.  */
default(parisizemax, 8*10^9);
read("/home/ubuntu/code/math-modular-sources/lattice/emn/central_ab.gp");
NN = #ACEN - 1;
acen = ACEN; bcen = BCEN;

/* normalised: ucen = 256^n b, vcen = 256^n a */
ucen = vector(NN+1, i, 256^(i-1)*bcen[i]);
vcen = vector(NN+1, i, 256^(i-1)*acen[i]);

fitrec(seqs, ord, dg, nmax) =
{
  my(rows = List(), unk = (ord+1)*(dg+1));
  for(s = 1, #seqs,
    my(W = seqs[s]);
    for(n = 0, nmax - ord,
      my(row = vector(unk));
      for(j = 0, ord, for(e = 0, dg, row[j*(dg+1)+e+1] = n^e * W[n+j+1]));
      listput(rows, row)
    )
  );
  my(mt = matconcat(Vec(rows)~));
  matker(mt);
}

print("=== fit U alone (256^n b_n = C(4n,2n)^2) ===");
for(ord = 1, 2, for(dg = 1, 6, my(k = fitrec([ucen], ord, dg, NN)); if(#k > 0, print("  order ",ord," deg ",dg,": nullspace dim ",#k); break(2))));

print("=== fit V alone (256^n a_n) ===");
found = 0;
for(ord = 1, 5, for(dg = 1, 12, my(k = fitrec([vcen], ord, dg, NN)); if(#k > 0, print("  order ",ord," deg ",dg,": nullspace dim ",#k); found=1; break(2))));
if(!found, print("  none with order<=5, deg<=12"));

print("=== fit the PAIR (U,V) jointly ===");
found = 0;
for(ord = 2, 5, for(dg = 1, 14, my(k = fitrec([ucen,vcen], ord, dg, NN)); if(#k > 0, print("  order ",ord," deg ",dg,": nullspace dim ",#k, "  kernel: ", k); found=1; break(2))));
if(!found, print("  NONE with order<=5, deg<=14"));
quit;
