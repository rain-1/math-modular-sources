default(parisizemax, 6000000000);
read("sc_rows.gp");
\\ formal test of a_{m p^s} = a_{m p^{s-1}} mod p^{c*s + d}
sup(nm, cf, r, N, spec) = {
  my(a = genseq(cf, r, [1], N));
  print("### ", nm, "  N=", N);
  for(i = 1, #spec,
    my(p = spec[i][1], c = spec[i][2], d = spec[i][3]);
    my(nt = 0, nf = 0, fb = 0);
    my(s = 1, ps = p);
    while(ps <= N,
      for(m = 1, N\ps,
        nt++;
        my(z = a[m*ps+1] - a[m*(ps\p)+1]);
        if(z % p^(c*s+d) != 0, nf++; if(fb==0, fb=[m,s,if(z==0,"0",valuation(z,p)), c*s+d]));
      );
      s++; ps *= p;
    );
    print("   p=", p, " modulus p^(", c, "s+", d, ") : ", nt, " tests, ", nf, " failures", if(nf, Str("  first (m,s,v_p,needed)=", fb), ""));
  );
};
sup("R1 Apery zeta(3) [Beukers-Coster control]", R1cf, 2, 3000, [[5,3,0],[7,3,0],[13,3,0],[2,3,-1],[3,3,-1]]);
sup("R2 Zagier C", R2cf, 2, 3000, [[2,2,0],[5,2,0],[7,2,0],[13,2,0],[3,2,0]]);
sup("R3 AZ eta", R3cf, 2, 3000, [[2,3,-2],[3,3,0],[7,3,0],[13,3,0],[5,3,0]]);
sup("R4 Sym^3 Zagier E", R4cf, 6, 3000, [[2,2,1],[2,2,0]]);
sup("R5 AESZ 207", R5cf, 4, 3000, [[3,2,0],[5,2,0],[7,2,0],[13,2,0],[2,3,0],[2,3,6]]);
