# herfurtner — scripts for `consolidation/HERFURTNER_CLASSIFICATION.md`

Classification of second-order Apery-like rows
`(n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1}` (P, Q quadratic) by the local
exponents of their Picard-Fuchs operator, matched against Herfurtner's list of
rational elliptic surfaces over P^1 with four singular fibres.

## Theory / tables

    01_exponents.py     exponent dictionary: rho_1, rho_2 at the finite singular
                        points, exponents at infinity, delta_infty, PSL_2(Z)
                        admissibility, Gauss-Bonnet signature.  Runs on the
                        known rows (Zagier's six, the nine Sym^1 square roots).
    04_herfurtner.py    Herfurtner Table 3 (the 38 rigid configurations with base
                        points, and the 18 one-parameter families), Euler-number
                        check, the cross-ratio invariant I = (1+z)^2/z = A^2/D,
                        and match_row().
    07_classtable.py    normalisation classes <-> configurations  -> out/classtable.txt

## The scan

    02_hscan.c          gcc -O2 -march=native -o 02_hscan 02_hscan.c
      ./02_hscan M J1 J2 AMIN AMAX BMAX CMAX N
        class (M; J1, J2):  Q(n) = C (Mn-J1)(Mn-J2),
                            P(n) = A(n^2 + (2M-J1-J2)/(2M) n) + B.
      Integrality tested division-free on U_n = u_n (n!)^2, prime by prime mod p^K.
    classes.txt         the 26 Kodaira-admissible classes
    run_scan.sh         drives the scan, 12 shards          (~75 min on 12 cores)
    run_a0.sh           the A = 0 pass

## Post-processing

    12_classsummary.py  per-class hit counts and family detection
    08_report.py        exact re-verification (n <= 300), companion b_n, sharp k,
                        lambda_i, score, cross-ratio match  -> out/rows_full.json
    09_table.py         markdown table of the non-degenerate rows
    13_fastfilter.py    all hits with |lambda_2| < 1 (arithmetic on (A,D) only)
    14_decay.py         deep analysis of those: k and score.  The only rows with
                        score > 0 and k >= 2 are Zagier D and Beukers'.

## The J-map test (is the projective monodromy inside PSL_2(Z)?)

    05_jtest.gp         builds the canonical nome q(t) from the Frobenius
                        recursion and tests whether J(gamma q^h) is a rational
                        function of t of degree <= 12.  A positive answer is a
                        certificate; it returns (h, deg J, gamma).
    05_run.gp 05_run2.gp 05_run3.gp   drivers for the known rows
    10_jtest_all.py     emits 11_jall.gp, the driver for every scanned row
                        (run:  gp -q 11_jall.gp < /dev/null > out/jall.log)

## Headline outputs

* `out/classtable.txt` — which Herfurtner configurations live in which
  normalisation class, with the finite list of admissible A^2/D per class.
* `out/jall.log` — the J-map verdicts.  Zagier's six (deg J = 12), the square
  roots of AZ(11,5,125) and AZ(9,3,-27) (deg J = 6), and two NEW rows
  (117,21,441) on I_1 I_7 II II (deg J = 8) and (72,6,108) on I_3 III III III
  (deg J = 3) are elliptic-surface rows; everything else is not.
