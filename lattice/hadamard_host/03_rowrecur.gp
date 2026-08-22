\\ 03_rowrecur.gp -- minimal recurrences for the two rows (in the aligned index n),
\\ and extension of the Nesterenko companion C_n from n<=98 to n<=NMAX.
default(parisizemax, 12000000000);
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read("/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/lib_fit.gp");
read(concat(outdir, "rows_raw.gp"));
NMAX = #AZ; NC = #BN0;
print("NMAX=", NMAX, "  NC=", NC);

\\ --- Zudilin 3-section: joint (Q_{3n}, P_{3n}) ---
fz = hh_fit("Zudilin 3-section (az,bz) = (Q_{3n},P_{3n})", [AZ, BZ], 6, 40, 240, NMAX);

\\ --- Nesterenko: an = 4 B_n alone (400 terms) ---
fn = hh_fit("Nesterenko a^N = 4 B_n", [AN], 6, 60, 300, NMAX);

\\ --- does the Nesterenko companion C_n satisfy the same recurrence? ---
{ if(type(fn) != "t_INT",
    my(r = fn[1], cs = fn[3], bad);
    bad = hh_verify([BN0], cs, r, 1, NC - r);
    print("");
    print("C_n satisfies the a^N recurrence on n=1..", NC-r, " : failures = ", bad);
    if(bad == 0,
      BN = hh_extend(BN0, cs, r, NMAX);
      write(concat(outdir, "rows_full.gp"), "BN = ", BN, ";");
      print("extended C_n to n=", NMAX, " and wrote out/rows_full.gp"),
      print("*** C_n does NOT satisfy it; fitting jointly instead");
      my(fj = hh_fit("Nesterenko joint (a^N,b^N)", [AN, BN0], 6, 60, NC, NC));
      if(type(fj) != "t_INT",
        BN = hh_extend(BN0, fj[3], fj[1], NMAX);
        write(concat(outdir, "rows_full.gp"), "BN = ", BN, ";");
        print("extended C_n to n=", NMAX, " (joint fit) and wrote out/rows_full.gp")))); }
\q
