# noncongruence_scan — scripts for `consolidation/NONCONGRUENCE_SCAN.md`

Two independent searches for Beukers-type irrationality rows
(`a_n = lambda^n [t^n] sqrt(F)`, `F` weight two with cuspidal divisor,
score `= log(1/|lambda_2|) - k`), plus the group-theoretic classification of the
possible host groups.

## Scan A — recurrence classes (ansatz-free, the load-bearing one)

    01_class_scan.c     gcc -O2 -o 01_class_scan 01_class_scan.c
      ./01_class_scan CLASS ALMIN ALMAX DEMAX GAMAX ZEMAX NTEST
        CLASS =  2,3,4,6   both finite singular points elliptic of order CLASS
                           (CLASS=2 is the R3 / free-integration / square-root class)
        CLASS =  0         both finite singular points are cusps (Zagier class)
        CLASS = -1         free common exponent difference: beta ranges freely and
                           epsilon = -2*delta*(alpha-beta)/alpha is derived
      Integrality is tested division-free on A_n = a_n (n!)^2, prime by prime mod p^K.
    04_score.py FILE OUT  Casoratian filter (Q(n) != 0 for n >= 1), measured k, score.
    11_summary.py         merge all *_scored.txt, dedupe, rank.
    08_period.py          Apery limit via the exact Casoratian series + identification.
    run_r3.sh run_e.sh run_free.sh run_more.sh run_more2.sh snake_wide.sh   drivers.

## Scan B — modular cross-check, eta-quotient pairs on Gamma_0(N), N <= 60

    02_eta_enum.py        Ligozat enumeration of (t, F) -> eta_pairs.json
    mkgp.py               eta_pairs.json -> lev/L<N>.gp
    mkjobs.py runjob.sh   sharding (jobs.txt); run with  xargs -P k -I{} ./runjob.sh "{}"
    03_modscan.gp         eta products; fitmin (exact) and fitminp (mod-p gated)
    03b_exact.gp          per level: lambda, minimal recurrence, char roots, k, score
    07_modtable.py        aggregation and ranking -> out/mod_table.json

## Host groups

    06_groups.py [nmax]   all genus-zero subgroups of PSL_2(Z) with exactly four
                          special points (index <= 12), with the Wohlfahrt congruence
                          test.  Needs sympy.

## The two classical Pade families (log at delta=16, arctan at delta=-27)

    13_ident234.gp        identification battery
    14_mirror234.gp       mirror map q = t exp(h/y1)  -- the modular/non-modular test
                          (integral for Apery's zeta(2) row and for Beukers' row;
                           NOT integral for either Pade family)
    15_verify234.gp       end-to-end Apery-criterion check on one row
    16_arctan.py          the arctan family: integrality (al = 18 mod 36), k=2,
                          and xi = -(1/(2 sqrt3)) arctan(6 sqrt3 / al) to 193 digits

## The level-5 Fricke row (the second non-congruence Beukers row)

    09_level5.gp          construction and verification (lambda = 2 minimal, k = 2 sharp)
    10_ident5.gp          xi to 260 digits + identification battery (negative)
    12_source5.gp         Theorem R3 source: Psi = f^3 t/2,  b_n = 2^n [t^n](f Theta)

All gp scripts are loaded as command-line arguments, e.g.

    gp -q 03_modscan.gp 03b_exact.gp lev/L6.gp driver.gp < /dev/null

(the `< /dev/null` matters: gp waits on stdin otherwise).
