# four_term — scripts for `consolidation/FOUR_TERM_SCAN.md`

Classification of **four-term** integral Apéry-like rows
`(n+1)^2 u_{n+1} = P(n)u_n - Q(n)u_{n-1} + R(n)u_{n-2}` (P, Q, R quadratic,
u_0 = 1, u_{-1} = u_{-2} = 0). The Picard-Fuchs operator

    L = theta^2 - t P(theta) + t^2 Q(theta+1) - t^3 R(theta+2)
      = t^2 Rc(t) D^2 + t Sc(t) D + t Vc(t),
      Rc = 1 - a t + d t^2 - g t^3

is **rank two with five singular points** 0, t_1, t_2, t_3, infinity -- the first
non-rigid case (two accessory parameters).  Sequel to `lattice/herfurtner/`.

## Theory / tables

    01_exponents.py   exponent dictionary, verified symbolically: the operator
                      identity, the D-form, T = Sc - t Rc', rho_i = -T(t_i)/(t_i Rc'(t_i)),
                      the Fuchs relation sum rho_i = s_1 + s_2 - 1, and the
                      equal-exponent conditions b = (1-rho)a, e = -2 rho d,
                      h = -(1+3 rho) g.
    02_fit.py         fits a four-term row to a sequence (9-parameter nullspace)
                      and runs the construction experiments: the gauge
                      y -> y/(1-nu t) and the signed binomial transform of a
                      three-term row, both of which produce four-term rows with an
                      APPARENT fifth singular point.
    04_classes.py     enumerates the Kodaira-admissible normalisation classes
                      (rho; M, j1, j2) with R(n) = C(Mn-j1)(Mn-j2).

## The scan

    03_fscan.c        gcc -O3 -march=native -o 03_fscan 03_fscan.c
      ./03_fscan RN RD M J1 J2 AMIN AMAX CMAX DMAX FMAX GMAX N [STRIDE SHIFT]
        rho = RN/RD;  P = a(n^2+(1-rho)n)+c,  Q = d(n^2-2 rho n)+f,
        R = C(Mn-J1)(Mn-J2);  loops a<=AMAX, |c|<=CMAX, |d|<=DMAX, |f|<=FMAX,
        |C|<=GMAX;  integrality prime-tested to n = N.
        DMAX = -1 and GMAX = -1 select the POSITIVE-SCORE WINDOW: |d| <= 2a+5,
        |C| <= (a+2)/M^2+1, which provably contains every row with |lambda_2| < 1.
        DMAX = -2 and GMAX = -2 select the POSITIVE-SCORE window: |d| <= 0.74(a+1)+2,
        |C| <= (0.136(a+1)+1)/M^2, which provably contains every row with
        |lambda_2| < 1/e, hence every row that could have score > 0 when k >= 1.
      Filters inside the scanner: cubic discriminant != 0; terminating sequences;
      constant-coefficient rows (P,Q,R all proportional to (n+1)^2); the exact
      analytic congruences 4|U_2, 36|U_3, 576|U_4 solved for C by CRT; a
      mask-arithmetic 2-adic prefilter to n=14 and a 64-bit prefilter for p<=11.
    classes.txt classes2.txt classes3.txt classes_all.txt classes_ord.txt  class lists
    run_scan.sh       ./run_scan.sh    TAG AMAX CMAX DMAX FMAX GMAX N NJOBS [classfile]
    run_window.sh     ./run_window.sh  TAG AMAX CMAX FMAX N NJOBS [classfile]   (|lambda_2|<1)
    run_pwindow.sh    ./run_pwindow.sh TAG AMAX CMAX FMAX N NJOBS [classfile]   (score>0)

## Post-processing

    05_report.py      exact re-verification to n=120; characteristic roots and
                      discriminant; the FROBENIUS OBSTRUCTION (apparent-singularity
                      test) at every t_i with rho in Z and at infinity when
                      delta_infty in Z; companion b_n, sharp k, Apery limit xi to
                      >= 60 digits with adaptive length, score log(1/|lambda_2|)-k.
                      `python3 05_report.py [full] < hits`
    08_analyse.py     deduplication: REPEATED / DISGUISED (apparent singularity =
                      four-point local system = gauge or cusp-move image of a
                      three-term row) / RESCALED (u_n -> mu^n u_n) / CANDIDATE.
                      `python3 08_analyse.py out/analysis.json out/*.txt`
    09_full.py        full invariants of the CANDIDATE rows (k, xi, score,
                      Moebius invariants I_1 = a^3/g, I_2 = a^2/d).

## The J-map test (deg J <= 18 by Riemann-Hurwitz, Theorem F3)

    06_jtest.gp       jtest5(rn,rd,M,j1,j2,a,c,d,f,C) -> [h, deg J, gamma, nfit, nextra]
                      or 0.  Canonical nome from the four-term Frobenius recursion
                      (u_n and g_n), then J(gamma q^h) fitted as U/V with
                      deg <= 18, filtered mod 2^61-1 and certified exactly over Q.
    10_jrun.py        emits the driver for every CANDIDATE row.
    13_jdetail.gp     jdetail(...) prints U, V, the factorisation of V (= the I_n
                      fibres and their n), the pole order at infinity, and the
                      factorisations of U and U - 1728 V (= the j=0 and j=1728
                      ramification) -- i.e. the complete Kodaira fibre configuration.

## Identification of the Apery limits

    07_ident.gp       reads out/xi.txt ("label xi") and runs lindep against
                      zeta(2), L(2,chi_D) for |D| <= 24, log^2, Pi*log, arctan,
                      Li_2, Gamma(1/4)^4/Pi, Gamma(1/3)^6/Pi^2 and products with
                      1/sqrt(m).
    12_polylog.gp     weight-two polylogarithm battery adapted to the row's own
                      singular points, with a greedy independent-subset reduction
                      of the basis (necessary: Li_2(1/2), Li_2(-1), Li_2(1/4)
                      satisfy classical relations that otherwise swamp lindep).

## Reproduction

    python3 01_exponents.py
    python3 02_fit.py
    python3 04_classes.py
    gcc -O3 -march=native -o 03_fscan 03_fscan.c
    ./run_scan.sh   G  200 40 1200 80 120 24 12 classes.txt    # general census
    ./run_window.sh W  400 60      120     24 12 classes.txt   # positive-score window
    python3 08_analyse.py out/analysis_all.json out/G*_*.txt out/P_*.txt out/W*_*.txt
    python3 09_full.py    out/analysis_all.json out/full.json
    python3 10_jrun.py    out/analysis_all.json 11_jall.gp
    gp -q -s 6000000000 11_jall.gp < /dev/null > out/jall.log
    gp -q -s 4000000000 07_ident.gp
