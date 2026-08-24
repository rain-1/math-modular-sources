# four_term_deep — scripts for `consolidation/FOUR_TERM_DEEP.md`

Deep extension of `lattice/four_term/` (the first four-term census).  Three new things:
the **mixed-exponent classes** (rho_1 != rho_2 = rho_3, Theorem D3), the **decayer
obstruction** (Theorems D1/D2), and a **five-term (six singular point) probe** (Theorem D5).

The headline results are in `consolidation/FOUR_TERM_DEEP.md`:
  * six genuine five-point four-term rows with Catalan / L-value Apery limits (Sec 6.4),
  * an integral six-point five-term row whose Apery limit is (sqrt2/3) L(g,2), the
    elliptic-K3 weight-3 newform value (Sec 7.3).

## Theory / checks

    01_mixcheck.py     Theorem D3 (mixed normalisation forms) against the Theorem F1
                       exponent dictionary, 15 classes x 6 random parameter points.
    04b_fit5.py        the five-term framework via the binomial transform of the K3 row.
    04c_brute5.py      helper brute force for the five-term scanner.
    04d_check5.py      independent exact brute force; validates 04_fscan5 exactly.
    09_prop.gp         Theorems D1/D2: the census rows, and an exhaustive numerical check
                       over 37088 integer cubics with |lambda_2| < 1 (0 reducible).

## Scanners (C)

    02_fscan.c   equal-rho, deep box.  Same interface as ../four_term/03_fscan, plus:
                 - the analytic-congruence chain extended by 14400 | U_5 (U_5 is still
                   LINEAR in C), which thins the innermost AP by 1-2 orders of magnitude;
                 - the modular inverses for the 36- and 576-congruences hoisted out of the
                   (d,f) loops (they depend only on the class and on (a,c));
                 - a modular (2 x 61-bit) triviality test, safe for large |d|.
                 Verified to return bit-identical hit lists to ../four_term/03_fscan.
    03_fmix.c    mixed-rho.  Class (rho_p, rho_r; M, J1, J2); loops over the rational
                 characteristic root r, then a, C, c, with f put on an AP by 4|U_2, 36|U_3.
                 usage: 03_fmix RPN RPD RRN RRD M J1 J2 RMIN RMAX AMAX CMAX DMAX FMAX GMAX N
                 out:   RPN RPD RRN RRD M J1 J2 r a c d f C
    04_fscan5.c  five-term (six singular points).  Class (rho; M, J1, J2), free parameters
                 a,c,d,f,g,j,C.  A trailing 1 as the 17th argument keeps repeated-root rows
                 (flagged "# REPEATED") so the census can be counted.
                 usage: 04_fscan5 RN RD M J1 J2 AMIN AMAX CMAX DMAX FMAX GMAX JMAX KMAX N
                                  [STRIDE SHIFT NODISC]

    mkjobs.py / run_queue.sh   a-block ordered job queues (low-a prefix always complete)
    classes_eq.txt classes_prod.txt classes_rest.txt classes_mix.txt classes_mix2.txt
    classes_5.txt              class lists

## Post-processing

    06_analyse_deep.py  equal-rho AND mixed-rho: exact integrality to n=120, characteristic
                        roots, the Frobenius obstruction, k, xi, score, t -> -t
                        canonicalisation, verdicts REPEATED / DISGUISED / RESCALED /
                        CANDIDATE.   python3 06_analyse_deep.py [--full] out.json files...
    05b_report5.py      the same for five-term rows.
    12_window.py        pulls the |lambda_2| < 1 rows out of a window scan and labels them.
    14_f5row.py         full invariants of one five-term row (integrality, k, xi).
    19_newrows.py       batch: analyse -> CANDIDATE -> "label value" file for the battery.

## Periods

    08_fold.py     fold (connection) constants of a general four-term row: builds Rc, Sc, Vc
                   from the class, continues (A,A') and (B,B') to each singular point, and
                   matches to the local Frobenius basis (logarithmic for rho in Z,
                   {1, s^rho} otherwise).  Self-test reproduces the K3 row's published
                   constants to 210 digits.   python3 08_fold.py selftest|census|row ...
    11_foldmix.py  the same for mixed-exponent rows (per-point rho, rotated-path fallback).
    15_verify_cat.py  from-scratch verification of the six Catalan rows.
    16_catconfirm.gp  their identifications at 300-digit working precision.

## Identification

    07_catalan.gp  the Catalan-focused battery: 145 constants, lindep at 60 digits with
                   every relation verified by reconstruction to 55 digits; T1 single, T2
                   pairs over a greedily independent 22-element sublist, T3 explicit
                   (1,xi,G,zeta(2)) vectors, T4 algdep to degree 6.  Self-tested.
                   gp -q -s 2000000000 07_catalan.gp                 (default input)
                   ./18_batt.sh /abs/path/limits.txt                 (any input)
    10_catsrc_*.py the Zagier E cusp-move family and Proposition D4.
    13_lvl16_*.gp / 13_lvl12_*.gp  the modular Catalan hosts (Sec 7.4).
    17_arith_*     the sigma_p / p-adic arithmetic of the six Catalan rows.

## Data

`out/` holds every scan output (`DG_*`, `DW_*`, `MIX_*`, `MX2_*`, `F5*`), the analysis
JSONs, `fold_existing.txt`, `cat_rows_xi.txt`, `cat_verify.log`, and the battery logs
`catalan_ident_*.log`.
